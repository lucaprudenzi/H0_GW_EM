bayestar-realize-coincs \
`# Write output to coinc.xml.` \
-o coinc_10000.xml \
`# Use the injections and noise PSDs that we generated.` \
inj_10000.xml --reference-psd psd.xml \
`# Specify which detectors are in science mode.` \
--detector H1 L1 V1 K1 \
`# Optionally, add Gaussian noise (rather than zero noise).` \
--measurement-error gaussian-noise \
`# Optionally, adjust the detection threshold: single-detector` \
`# SNR, network SNR, and minimum number of detectors above` \
`# threshold to form a coincidence.` \
--snr-threshold 4.0 \
--net-snr-threshold 8.0 \
--min-triggers 2 \
--max-distance 400
`# Optionally, save triggers that were below the single-detector` \
`# threshold.` \
--keep-subthreshold
-j 64
