%% make correlation smoothing matrices 
dirname = '..\MATLAB\Data\saitamav2\';
load(strcat(dirname, 'ID.mat'));
load(strcat(dirname, 'data.mat'));
in = matfile(strcat(dirname, 'in.mat'));
param = matfile(strcat(dirname, 'param.mat'), 'Writable', true);

samples = unique([data.Sample]);
[spectra, spectrumsIdx, ~]=  unique(strcat({ID.Csvid}, {ID.T}));

wavelengthN = 81;
if (wavelengthN == 81)
    wavelength = 380:5:780;
    ver = '';
elseif (wavelengthN == 401)
    wavelength = 380:780;
    ver = '_401';
else 
    warning('Not supported length')
    return;
end 

%% case 'known correlation fixed cancer sample' % average xcorr of all fixed cancer measured data for fixed cancer sample

for i = 1:length(samples)
    %find the relative samples
    name = samples{i};
    spectrum = spectra( contains(spectra, name));
    spectrumIdx = spectrumsIdx( contains(spectra, name));
    
    r = zeros(81, length(spectrum));
    rnc = [];
    rnf = [];
    rnu = [];
    rcc = [];
    rcf = [];
    rcu = [];
    
    for j=1:length(spectrum)
        id = ID(spectrumIdx(j));
        r(:,j) = interp1(380:780, measuredSpectrumStruct(j).Spectrum, wavelength, 'nearest');
        if (id.IsNormal)
            if (id.IsCut)
                rnc = [rnc, r(:,j)];
            elseif (id.IsFixed)
                rnf = [rnf, r(:,j)];
            else 
                rnu = [rnu, r(:,j)];                
            end
        else 
            if (id.IsCut)
                rcc = [rcc, r(:,j)];
            elseif (id.IsFixed)
                rcf = [rcf, r(:,j)];
            else 
                rcu = [rcu, r(:,j)];                
            end
        end
    end
    
    
    param.(matlab.lang.makeValidName(['Mnc_', name])) = 1/size(rnc,2) * (rnc * rnc');
    param.(matlab.lang.makeValidName(['Mnf_', name])) = 1/size(rnf,2) * (rnf * rnf');
    param.(matlab.lang.makeValidName(['Mnu_', name])) = 1/size(rnu,2) * (rnu * rnu');

    param.(matlab.lang.makeValidName(['Mcc_', name])) = 1/size(rcc,2) * (rcc * rcc');
    param.(matlab.lang.makeValidName(['Mcf_', name])) = 1/size(rcf,2) * (rcf * rcf');
    param.(matlab.lang.makeValidName(['Mcu_', name])) = 1/size(rcu,2) * (rcu * rcu');
    
end


%% case 'known correlation cancer sample'
for i = 1:length(samples)
    % find the relative samples
    name = samples{i};
    spectrum = spectra( contains(spectra, name));
    spectrumIdx = spectrumsIdx( contains(spectra, name));
    
    r = zeros(81, length(spectrum));
    rn = [];
    rc = [];
    for j=1:length(spectrum)
        id = ID(spectrumIdx(j));
        r(:,j) = interp1(380:780, measuredSpectrumStruct(j).Spectrum, wavelength, 'nearest');
        if (id.IsNormal)
            rn = [rn, r(:,j)];
        else 
            rc = [rc, r(:,j)];
        end
    end
    
    % correlation of all the measured spectra for this sample
    param.(matlab.lang.makeValidName(['M_', name])) = 1/length(spectrum) * (r * r');
    
    % correlation of all cancerous spectra for this sample
    param.(matlab.lang.makeValidName(['Mc_', name])) = 1/size(rc,2) * (rc * rc');

    % correlation of all normal spectra for this sample
    param.(matlab.lang.makeValidName(['Mn_', name])) = 1/size(rn,2) * (rn * rn');
    
end

            
% %% case 'markovian'
%     rho = 0.98;
%     M_row = ones(1,wavelengthN);
%     for i=2:wavelengthN 
%         M_row(i) = M_row(i-1) * rho;
%     end
% 
%     M = zeros(wavelengthN); %first order Markov process covariance matrix
%     select = eye(1, wavelengthN);
%     D = zeros(wavelengthN,wavelengthN); % Q x lambdaN
%     for i=1:wavelengthN
%         M(i,:) = [ fliplr(M_row(1: i)) , M_row(2:wavelengthN-i+1)];
%         D(i,:) = circshift( select', (i-1)*wavelengthN); %diagonal with 1 in the position of components
%     end
%     param.(matlab.lang.makeValidName(['M_markov_98', ver])) = toeplitz(r);

%% case 'known correlation all data'
   % m = matfile('refestparam.mat');
    [meas, ~, idx, ~] = subset('measured', name,  'unique');
    avgMeasured = interp1(380:780, mean(meas), wavelength, 'nearest');
    z = xcorr(avgMeasured, avgMeasured, 'coeff');
    %z=xcorr(m.avgMeasured, m.avgMeasured, 'coeff');
    r = z(wavelengthN:end);
    param.(matlab.lang.makeValidName('M_measured_xcorr')) = toeplitz(r);

% %% case 'known correlation macbeth'
%     m = matfile('refestparam.mat');
%     z=xcorr(m.avgMeasuredMacbeth, m.avgMeasuredMacbeth, 'coeff');
%     r = z(wavelengthN:end);
%     param.(matlab.lang.makeValidName(['M_measured_xcorr_macbeth',ver])) = toeplitz(r);

%% case 'known correlation single'
%     %old implementation 
%     %     [~, Ks] = corrmtx(measuredReflectance, wavelengthN-1);
%     z=xcorr(spectrum, spectrum, 'coeff');
%     r = z(wavelengthN:end);
%     M =toeplitz(r);

%%     case 'KCor all same malignancy'
rn = [];
rc = [];
for i = 1:length(ID)
        ref = interp1(380:780, measuredSpectrumStruct(i).Spectrum, wavelength, 'nearest')';
        if (ID(i).IsNormal)
            rn = [rn, ref];
        else 
            rc = [rc, ref];
        end
end   
% correlation of all cancerous spectra
param.(matlab.lang.makeValidName('Mc')) = 1/size(rc,2) * (rc * rc');

% correlation of all normal spectra 
param.(matlab.lang.makeValidName('Mn')) = 1/size(rn,2) * (rn * rn');
    
%%     case 'KCor all same fixation'
rf = [];
ru = [];
rc = [];
for i = 1:length(ID)
        ref = interp1(380:780, measuredSpectrumStruct(i).Spectrum, wavelength, 'nearest')';
        if (ID(i).IsCut)
            rc = [rc, ref];
        elseif (ID(i).IsFixed)
            rf = [rf, ref];
        else 
            ru = [ru, ref];
        end
end   
% correlation of all fixed spectra
param.(matlab.lang.makeValidName('M_f')) = 1/size(rf,2) * (rf * rf');

% correlation of all unfixed spectra
param.(matlab.lang.makeValidName('M_u')) = 1/size(ru,2) * (ru * ru');

% correlation of all cut spectra
param.(matlab.lang.makeValidName('M_c')) = 1/size(rc,2) * (rc * rc');
  

    