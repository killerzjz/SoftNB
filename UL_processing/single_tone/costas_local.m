function[output] = costas_local(nb_input)
previous_sample = 0+0j;
phase = 0;
DDSPreviousInp = 0;
loopFilState = 0;
integFiltState = 0;
IntegratorGain = 1.3706e-06;
proportionalGain = 1.0278e-04;
err = [];
for i = 1:length(nb_input)
    %Find Phase Error(Phase Detector)
    phErr = sign(real(previous_sample)).*imag(previous_sample)...
              - sign(imag(previous_sample)).*real(previous_sample);
    err = [err phErr];
    output(i) = nb_input(i)*exp(1j*phase);

    %Loop Filter
    loopFilout = loopFilState+IntegratorGain*phErr;
    loopFilState = loopFilout;

    % Direct digital synthesizer implemented as an integrator
    DDSOut = DDSPreviousInp + integFiltState;
    integFiltState = DDSOut;
    DDSPreviousInp = phErr*proportionalGain+loopFilout;

    phase = -1 * DDSOut;
    previous_sample = output(i);
end
end