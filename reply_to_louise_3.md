Subject: Re: 2nd-gen biofuel rampdown after 2450 (VL)

Hi Louise,

Good catch. It's **VL-only** — H and M don't show it — and it comes from the
input emissions, not our extension code.

We scale the biofuel fractions (crpbiof and the 2nd-gen crpbf_c3per/c4per) by a
BECCS factor = each scenario's BECCS CO₂ relative to its 2100 value, taken
straight from the harmonised emissions CSV. For H and M that factor is flat at
1.0 (their BECCS is ~zero / constant), so the biofuel fields are constant after
2100. VL is the only scenario with a time-varying BECCS trajectory, and in the
harmonised timeseries VL's BECCS holds near −2760 Mt CO₂/yr to ~2450 and then
ramps down to essentially zero by 2500:

    2450: −2758   2470: −958   2490: −97   2500: −1.0 (Mt CO₂/yr)

So the biofuel fractions faithfully follow that end-of-horizon BECCS phase-out.
It shows up most clearly in the 2nd-gen crpbf fields because they're small and
never hit the 1.0 cap, so they track the factor directly. Plot attached (H sits
under M at 1.0).

Bottom line: it's a genuine feature of the REMIND-MAgPIE VL emissions extension
that we propagate correctly — not a bug, and it doesn't affect H or M. The only
judgement call is whether that 2450–2500 BECCS drawdown is intended in the
harmonised emissions or something we'd rather cap/hold flat; happy to do the
latter as a modelling choice if you'd prefer.

Best,
Ben

---
Attachment: output/beccs_scaling_diagnostic.png
