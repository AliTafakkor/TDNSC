function modulationEnergy = modulation_energy(coch, P, ti, envelopes)

if nargin < 4
    envelopes = false;
end

spec_mod_rates = P.spec_mod_rates;
pos_temp_mod_rates = P.temp_mod_rates(P.temp_mod_rates>0);
temp_mod_rates = pos_temp_mod_rates;
ti = 1:size(coch,1);

% 2D fourier transform of the padded cochleogram
FT_padded_coch = fft2(pad_coch(coch,P));

% dimensions of frequency-padded cochleogram
[T,F] = size(FT_padded_coch);

n_spec_mod_rates = length(spec_mod_rates);
n_temp_mod_rates = length(temp_mod_rates);
modulationEnergy = nan(n_spec_mod_rates, n_temp_mod_rates);
for i = 1:n_spec_mod_rates
    for j = 1:n_temp_mod_rates
        
        % transfer function of spectrotemporal filter
        Hts = filt_spectemp_mod(...
            spec_mod_rates(i), temp_mod_rates(j), F, T, P,...
            [], [], [], [], [], [], [], ...
            P.spec_BW, P.temp_BW, P.spec_wavelet, P.temp_wavelet, ...
            P.spec_random_phase, P.temp_random_phase, ...
            P.spec_random_filt, P.temp_random_filt, P.random_seed);
        
        % whether or not to compute envelopes of the signal
        if envelopes
            
            if ~isnan(temp_mod_rates(j)) % compute temporal envelopes if possible
                Hts = analytic_from_spectrum_2D(Hts, 1);
            elseif ~isnan(spec_mod_rates(i)) % otherwise compute spectral envelopes
                Hts = analytic_from_spectrum_2D(Hts, 2);
            else
                error('Envelopes cannot be computed for filter without any modulation');
            end
            
            % apply transfer function, absolute value of analytic signal
            filtcoch_padded = abs(ifft2(FT_padded_coch .* Hts));
            filtcoch = remove_pad(filtcoch_padded, P);
            
        else
            
            % apply transfer function
            filtcoch_padded = real(ifft2(FT_padded_coch .* Hts));
            filtcoch = remove_pad(filtcoch_padded, P);
        
        end

        modulationEnergy(i,j) = sum(filtcoch(ti,:).^2, 'all');
        
    end
end
