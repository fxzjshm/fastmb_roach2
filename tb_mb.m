%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Validate Single Beam Model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sys = 'mb1k';
VERBOSE = 1;
fprintf('\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Spectrometer
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% collect model configuration parameters
pfbsize = eval(get_param([sys, '/u0/x4_spec'], 'FFTSize'));
ntaps = eval(get_param([sys, '/u0/x4_spec/pfb_fir_real'], 'TotalTaps'));
pfb_wintype = get_param([sys, '/u0/x4_spec/pfb_fir_real'], 'WindowType');
spec_bp = eval(get_param([sys, '/u0/x4_vacc/vacc_ufix'], 'out_bin_pt'));
fftsize = pfbsize;
nfft = 2^fftsize;
nbins = nfft / 2;

% collect simulation parameters
gain = eval(get_param([sys, '/u0/sim_gain'], 'value'));
gain0 = bitand(gain, 2^16-1, 'uint32');
gain1 = bitshift(gain, -16, 'uint32');
fftshift = eval(get_param([sys, '/u0/sim_fftshift'], 'value'));
acclen = eval(get_param([sys, '/u0/sim_acclen'], 'value'));
bsval = eval(get_param([sys, '/u0/sim_bitsel'], 'value'));
% bitsel is bit pattern (0-1, 2-3, 4-5, 6-7) map to 4 stokes
bitsel = zeros(4, 1);
for i = 1:4
    bitsel(i) = bitget(bsval, i*2) * 2 + bitget(bsval, i*2-1);
end

in_scope = ScopeIn;
fft_scope = ScopeFFT;
acc_scope = ScopeACC;
accf_scope = ScopeACCF;
pkt0_scope = ScopePkt0;
pkt1_scope = ScopePkt1;

% read in adc signals
first = find(in_scope.signals(1).values, 1, 'last') + 1;
total = length(in_scope.signals(1).values);
nelems = total - first + 1;
all0 = [in_scope.signals(2:5).values]';
in0 = reshape(all0(:, first:end), nelems * 4, 1);
all0 = [in_scope.signals(6:9).values]';
in1 = reshape(all0(:, first:end), nelems * 4, 1);
fprintf('Read in %d adc signals starting at %d\n', nelems * 4, first);

% read in powered fft
first = find(fft_scope.signals(1).values, 1, 'last') + 1;
nelems = length(fft_scope.signals(1).values) - first + 1;
fftout = zeros(nelems*2, 4);
for p = 1:4
    all0 = [fft_scope.signals((p-1)*2+2:p*2+1).values]';
    fftout(:, p) = reshape(all0(:, first:end), nelems * 2, 1);
end
fprintf('Read in %d fft spectrum starting at %d\n', nelems * 2, first);

nspecs = floor(nelems*2 / nbins);

% split fft scope
stokes_s = cell(4,1);
for p = 1:4, stokes_s{p} = zeros(nbins, nspecs); end
for i = 1:nspecs
    for p = 1:4
        stokes_s{p}(:, i) = fftout((i-1)*nbins+1:i*nbins, p);
    end
end

% read in vacc
first = first + 8;  % vacc has 6+1+1 clocks delay
valid = find(acc_scope.signals(1).values(first:end)) + (first - 1);
nelems = length(valid);
nacc = floor(nelems*2 / nbins);
acc_s = cell(4,1);
for p = 1:4
    lo = (p-1)*2 + 1 + 1;
    hi = p*2 + 1;
    all0 = [acc_scope.signals(lo:hi).values];
    acc = transpose(all0(valid, :));
    acc_s{p} = reshape(acc(:, 1:nbins*nacc/2), nbins, nacc);
end

% full vacc
accf_s = cell(4,1);
for p = 1:4
    lo = (p-1)*2 + 1 + 1;
    hi = p*2 + 1;
    all0 = [accf_scope.signals(lo:hi).values];
    acc = transpose(all0(valid, :));
    accf_s{p} = reshape(acc(:, 1:nbins*nacc/2), nbins, nacc);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Calculate spectrum from ScopeIn
%
pfb = polyfb(nfft, ntaps, pfb_wintype);
pfbout0 = pfb.apply(in0(1:nfft*(ntaps+nspecs-1)));
pfbout1 = pfb.apply(in1(1:nfft*(ntaps+nspecs-1)));
pfbout0 = reshape(pfbout0, nfft, nspecs);
pfbout1 = reshape(pfbout1, nfft, nspecs);

% calculate number of fft shift
value = fftshift;
shiftcount = 0;
for i = 1 : fftsize
    if bitand(value, 1) ~= 0
        shiftcount = shiftcount + 1;
    end
    value = bitshift(value, -1);
end

stokes_c = cell(4,1);
for p = 1:4, stokes_c{p} = zeros(nbins, nspecs); end
% don't know where the extra divided by 2 comes from
scale0 = (gain0 / 2^8) / 2^shiftcount / 2;
scale1 = (gain1 / 2^8) / 2^shiftcount / 2;
for i = 1:nspecs
    fftout0 = fft(pfbout0(:, i)) * scale0;
    fftout1 = fft(pfbout1(:, i)) * scale1;
    stokes_c{1}(:, i) = abs(fftout0(1:nbins)) .^ 2;
    stokes_c{2}(:, i) = abs(fftout1(1:nbins)) .^ 2;
    fftout_cross = fftout0(1:nbins) .* conj(fftout1(1:nbins));
    stokes_c{3}(:, i) = real(fftout_cross);
    stokes_c{4}(:, i) = imag(fftout_cross);
end

for p = 1:4
    rms_s = rms(stokes_s{p});
    max_s = max(stokes_s{p});
    rms_c = rms(stokes_c{p});
    max_c = max(stokes_c{p});
    delta = abs(stokes_c{p} - stokes_s{p});
    rms_delta = rms(delta);
    max_delta = max(delta);
    rel_delta = rms_delta ./ rms_c;

    fprintf('\nFFT spectrum statistics of polarization %d:\n', p);
    if VERBOSE
        fprintf('\n');
        fprintf('  Nr   rms_s     rms_c   rms_diff    max_s     max_c   max_diff  rel_diff\n');
        fprintf('---------------------------------------------------------------------------\n');
        for i = 1 : nspecs
            fprintf('%3d  %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e\n', ...
                     i, rms_s(i), rms_c(i), rms_delta(i), ...
                     max_s(i), max_c(i), max_delta(i), rel_delta(i));
        end
        fprintf('\n');
    end
    fprintf('Max relative error %g, rms relative error %g\n', ...
            max(rel_delta), rms(rel_delta));
end

acc_c = cell(4,1);
accf_c = cell(4,1);
stride = acclen + 1;
for p = 1:4
%     scale = 2^(spec_bp - 8*bitsel(p));
    for i = 1:nacc
        acc = sum(stokes_c{p}(:, (i-1)*stride+1 : i*stride), 2);
        accf_c{p}(:, i) = acc * 2^spec_bp;
        acc_c{p}(:, i) = accf_c{p}(:, i) / 2^(8*bitsel(p));
    end
end

for p = 1:4
    rms_s = rms(accf_s{p});
    max_s = max(accf_s{p});
    rms_c = rms(accf_c{p});
    max_c = max(accf_c{p});
    delta = abs(accf_c{p} - accf_s{p});
    rms_delta = rms(delta);
    max_delta = max(delta);
    rel_delta = rms_delta ./ rms_c;

    fprintf('\nACC FULL spectrum statistics of polarization %d:\n', p);
    if VERBOSE
        fprintf('\n');
        fprintf('  Nr   rms_s     rms_c   rms_diff    max_s     max_c   max_diff  rel_diff\n');
        fprintf('---------------------------------------------------------------------------\n');
        for i = 1 : nacc
            fprintf('%3d  %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e\n', ...
                     i, rms_s(i), rms_c(i), rms_delta(i), ...
                     max_s(i), max_c(i), max_delta(i), rel_delta(i));
        end
        fprintf('\n');
    end
    fprintf('Max relative error %g, rms relative error %g\n', ...
            max(rel_delta), rms(rel_delta));
end


for p = 1:4
    rms_s = rms(acc_s{p});
    max_s = max(acc_s{p});
    rms_c = rms(acc_c{p});
    max_c = max(acc_c{p});
    delta = abs(acc_c{p} - acc_s{p});
    rms_delta = rms(delta);
    max_delta = max(delta);
    rel_delta = rms_delta ./ rms_c;

    fprintf('\nACC QUANT spectrum statistics of polarization %d:\n', p);
    if VERBOSE
        fprintf('\n');
        fprintf('  Nr   rms_s     rms_c   rms_diff    max_s     max_c   max_diff  rel_diff\n');
        fprintf('---------------------------------------------------------------------------\n');
        for i = 1 : nacc
            fprintf('%3d  %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e\n', ...
                     i, rms_s(i), rms_c(i), rms_delta(i), ...
                     max_s(i), max_c(i), max_delta(i), rel_delta(i));
        end
        fprintf('\n');
    end
    fprintf('Max relative error %g, rms relative error %g\n', ...
            max(rel_delta), rms(rel_delta));
end

% check quantization
quant_err = cell(4, 1);
for p = 1:4
    quant_err{p} = floor(accf_s{p} ./ 2^(8*bitsel(p))) - acc_s{p};
    for i = 1:nacc
        nz = find(quant_err{p}(:,i));
        if ~isempty(nz)
            warning('Error found in bit selection for pol %d acc %d:', p, i);
            disp(nz);
        end
    end
end

% read in spectrometer packets
first = find(pkt0_scope.signals(1).values, 1, 'last') + 1;
valid = find(pkt0_scope.signals(10).values(first:end)) + (first - 1);
eof = find(pkt0_scope.signals(11).values(first:end)) + (first - 1);
all0 = [pkt0_scope.signals(2:9).values];
all1 = [pkt1_scope.signals(2:9).values];
payload0 = [];
payload1 = [];
pkt_begin = 1;
for i = 1:length(eof)
    pkt_end = find(valid == eof(i));
    sn0 = all0(valid(pkt_begin), :);
    sn1 = all1(valid(pkt_begin), :);
    if sn0(1)+1 ~= i || sn1(1)+1 ~= i
        error('Wrong serial number!')
    end
    payload0 = [payload0 all0(valid(pkt_begin+1:pkt_end), :)'];
    payload1 = [payload1 all1(valid(pkt_begin+1:pkt_end), :)'];
    pkt_begin = pkt_end + 1;
end

% convert to fixed value
ufix8_to_fix8 = @(x) (1-bitget(x,8)) .* x + -bitget(x,8) .* (bitcmp(x,'uint8') + 1);
% payload0 = ufix8_to_fix8(payload0);
payload1 = ufix8_to_fix8(payload1);

nspecs = floor(length(payload0) * 4 / nbins);    % 4 parallel
vlen = nspecs*nbins/4;
acc_p = cell(4,1);
% aa = vertcat(payload0(1:2, 1:vlen), payload0(5:6, 1:vlen));
% acc_p{1} = reshape(aa, nbins, nspecs);
% bb = vertcat(payload0(3:4, 1:vlen), payload0(7:8, 1:vlen));
% acc_p{2} = reshape(bb, nbins, nspecs);
% cr = vertcat(payload1(1:2, 1:vlen), payload1(5:6, 1:vlen));
% acc_p{3} = reshape(cr, nbins, nspecs);
% ci = vertcat(payload1(3:4, 1:vlen), payload1(7:8, 1:vlen));
% acc_p{4} = reshape(ci, nbins, nspecs);
aa = payload0(1:2:8, 1:vlen);
acc_p{1} = reshape(aa, nbins, nspecs);
bb = payload0(2:2:8, 1:vlen);
acc_p{2} = reshape(bb, nbins, nspecs);
cr = payload1(1:2:8, 1:vlen);
acc_p{3} = reshape(cr, nbins, nspecs);
ci = payload1(2:2:8, 1:vlen);
acc_p{4} = reshape(ci, nbins, nspecs);


% check packetizer
pkt_err = cell(4, 1);
for p = 1:4
    pkt_err{p} = acc_p{p} - acc_s{p};
    for i = 1:nacc
        nz = find(pkt_err{p}(:,i));
        if ~isempty(nz)
            warning('Error found in packetization for pol %d acc %d:', p, i);
            disp(nz);
        end
    end
end

for p = 1:4
    rms_s = rms(acc_p{p});
    max_s = max(acc_p{p});
    rms_c = rms(acc_c{p}(:,1:nspecs));
    max_c = max(acc_c{p}(:,1:nspecs));
    delta = abs(acc_c{p}(:,1:nspecs) - acc_p{p});
    rms_delta = rms(delta);
    max_delta = max(delta);
    rel_delta = rms_delta ./ rms_c;

    fprintf('\nPKT spectrum statistics of polarization %d:\n', p);
    if VERBOSE
        fprintf('\n');
        fprintf('  Nr   rms_s     rms_c   rms_diff    max_s     max_c   max_diff  rel_diff\n');
        fprintf('---------------------------------------------------------------------------\n');
        for i = 1 : nspecs
            fprintf('%3d  %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e\n', ...
                     i, rms_s(i), rms_c(i), rms_delta(i), ...
                     max_s(i), max_c(i), max_delta(i), rel_delta(i));
        end
        fprintf('\n');
    end
    fprintf('Max relative error %g, rms relative error %g\n', ...
            max(rel_delta), rms(rel_delta));
end
