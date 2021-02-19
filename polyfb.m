classdef polyfb
    
    properties (SetAccess = private)
        nfft
        ntaps
        coeff
    end

    methods
        function obj = polyfb(nfft, ntaps, window_type)

            if nfft <= 0 || log2(nfft) ~= nextpow2(nfft)
                err = MException('PolyFB:OutOfRange', ...
                        'nfft must be positive integer power of 2');
                throw(err);
            end
            obj.nfft  = nfft;

            if ntaps <= 0 || mod(ntaps, 1) ~= 0
                err = MException('PolyFB:OutOfRange', ...
                        'nfft must be positive integer');
                throw(err);
            end
            obj.ntaps = ntaps;

            coeff_arr = sinc(-ntaps/2 : 1/nfft : ntaps/2-1/nfft)' ...
                            .* window(window_type, nfft*ntaps);
            obj.coeff = reshape(coeff_arr, nfft, ntaps);
        end
        
        function res = apply(obj, arr)
            nblk = length(arr)/obj.nfft - obj.ntaps + 1;
            res = zeros(obj.nfft * nblk, 1);
            blk = reshape(arr, obj.nfft, nblk + obj.ntaps - 1);
            for n = 1 : nblk
                for i = 1 : obj.nfft
                    res((n-1)*obj.nfft + i) = blk(i, n : n+obj.ntaps-1) * obj.coeff(i, :)';
                end
            end
        end
    end

end
