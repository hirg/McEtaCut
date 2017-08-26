# McEtaCut

### McPhiDecay
> `McPhiEta.C`
> - generate phi-meson through MC with AMPT eta, STAR published spectra and v2 
> - let phi-meson decay with PHYTIA and boost K+ back to phi-meson rest frame
> - correlate phi-meson with fixed event plane (Psi = 0)
> - extract rho_00 vs. eta_max

> `McPhiEtaBoost.C`
> - generate phi-meson at rest through MC
> - let phi-meson decay with PHYTIA
> - correlate phi-meson with fixed event plane (Psi = 0)
> - extract rho_00 vs. eta_max and compare with analytic calculation

### McLambdaDecay
> `McLambdaEta.C`
> - generate Lambda/anti-Lambda through MC with AMPT eta, STAR published spectra and v2 
> - let Lambda/anti-Lambda decay with PHYTIA and boost proton back to Lambda rest frame
> - correlate Lambda with fixed event plane (Psi = 0)
> - extract pH vs. eta_max
> - pH can be set to any value
> - use submit/McLambdaEta.sh to submit jobs to cumulate more statistics

> `McLambdaEtaBoost.C`
> - generate Lambda/anti-Lambda at rest through MC
> - let Lambda/anti-Lambda decay with PHYTIA
> - correlate Lambda with fixed event plane (Psi = 0)
> - extract pH vs. eta_max and compare with analytic calculation

### Utility
> - constant used in `McPhiDecay` and `McLambdaDecay`
> - functions used in `McPhiDecay` and `McLambdaDecay`

### figures
> - place to save QA and final plots
