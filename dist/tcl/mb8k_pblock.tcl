
proc new_pblock_regexp {pblk args} {
    if {[get_pblocks -q $pblk] != ""} { delete_pblock $pblk }
    create_pblock $pblk
    foreach pattern $args {
        add_cells_to_pblock $pblk [get_cells -hier -regexp $pattern]
    }
    place_pblocks -util 75 $pblk
    set_property gridtypes {RAMB36 RAMB18 DSP48 SLICE} [get_pblocks $pblk]
}

foreach zdok {zdok0 zdok1} {
	new_pblock_regexp ${zdok} .+?/${zdok}.+?/scope.+ .+?/${zdok}.+?/x4_tvg.+
}

foreach u {u0 u1} {
	new_pblock_regexp ${u}_pfbfir .+?/${u}.+?/x4_spec.+?/pfb_fir_real.+
	new_pblock_regexp ${u}_fft_biplex .+?/${u}.+?/x4_spec.+?/fft_wideband_real.+?/fft_biplex_real_4x.+
	new_pblock_regexp ${u}_fft_direct .+?/${u}.+?/x4_spec.+?/fft_wideband_real.+?/fft_direct.+
	new_pblock_regexp ${u}_fft_unscrambler .+?/${u}.+?/x4_spec.+?/fft_wideband_real.+?/fft_unscrambler.+
	new_pblock_regexp ${u}_other .+?/${u}.+?/x4_spec.+?/prescale.+ .+?/${u}.+?/x4_spec.+?/stokes.+
	new_pblock_regexp ${u}_vacc .+?/${u}.+?/x4_vacc.+
}

