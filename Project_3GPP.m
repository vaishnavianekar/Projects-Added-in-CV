fc = 3.5e9;
c = 3e8;
BW = 180e3; % 1 PRBs NB in hz
Pt_dbm = -5:5:35; %tx power 
Gt = 0; Gr =0; % antenna gains
NF_db = 7; %noise factor
kT_dbmhzz = -174;
nframes = 10e5; %no of frames for monte carlo simulation
Nsysm = 1024; %no of symbols in each frame
dmin = 10; %min distance in m
dmax = 100; % max distance in m

%lets define the pathloss model for 1m CI model
PL_1m = 20*log(4*pi*fc/c); %in dB

%pathloss exponenets and log normal shadowing std
n_LOS = 2 ; 
n_NLOS = 3.75;
sigma_LOS = 3;
sigma_NLOS = 8;

% rician parameter for LOS
K_db_LOS = 8;
K_LOS = 10^(K_db_LOS/10);

%Thresholds for AMC scheme
th_qpsk_to_16qam = 8 ; %if snr>= 8db switch to 16qam
th_16qam_to64qam = 18; %if snr >= 18db swicth to 64qam
%for lowest MCS use QPSK
fixedmodulation = 4; %M
modNames= containers.Map([1,2,4,16],['BPSK', 'pi/2-BPSK','QPSK',"16QAM"]); %mappnig numbers to names

%bits per symbol
bps = @(M) log2(M); %one line function in matlab bps is the func here

% variables to store the results
avgRxSNR_db = zeros(size(Pt_dbm));
ber_fixed = zeros(size(Pt_dbm));
ber_amc = zeros(size(Pt_dbm));
thorughput_fixed = zeros(size(Pt_dbm)); %bits/sec (nsymbs*bps / frametime)
throughput_amc = zeros(size(Pt_dbm));

%LOS probbaility function
pLOS_UMi = @(d) (d<= 18).*1+ (d>18).*(18./d + exp(-d/36).*(1-18./d));

%for modulation we will use built in functions
%main loop
for ip = 1:length(Pt_dbm)
    pt_dbm = Pt_dbm(ip);
    biterrors_fixed =0;
    biterrors_amc = 0;
    totalbits_fixed = 0;
    totalbits_amc =0;
    snr_acc =0;

    for frame = 1:nframes
        %random distance
        d = dmin + (dmax-dmin)*rand();
        %LOS prob choosing and selction
        pLOS = pLOS_UMi(d);
        isLOS = rand() < pLOS;

        %now choosing the pathloss exponent model+ shadowing

        if isLOS
            n = n_LOS;
            sigma = sigma_LOS;
        else
            n = n_NLOS;
            sigma = sigma_NLOS;
        end
        % setting the pathloss
        PLdb = PL_1m + 10*n*log10(max(d,1)) + sigma*randn();
        %received power 
        Pr_dbm = pt_dbm + Gt- Gr - PLdb;
        %linear powers
        Pr_w = 10.^((Pr_dbm -30)/10);
        noiseP_dbm = kT_dbmhzz + 10*log10(BW)+NF_db;
        N0_W = 10.^((noiseP_dbm-30)/10);

        %instantaneous SNR
        SNR_lin =  Pr_w/N0_W;
        SNR_db = 10*log10(SNR_lin);
        snr_acc = snr_acc + SNR_db;

        %fixed modulation part
        M_fixed = fixedmodulation;
        k_fixed = bps(M_fixed);
        % AMC: pick by threshold
        if SNR_db < th_qpsk_to_16qam
            M_amc = 4; %QPSK
        elseif SNR_db < th_16qam_to64qam
            M_amc = 16; % 16-QAM
        else
            M_amc = 64; % 64-QAM;
        end
        k_amc = bps(M_amc);

        %fixed modulation part
        Nsym = Nsysm;
        bits_fixed = randi([0 1], Nsym*k_fixed,1);
        sym_fixed = qammod(bits_fixed,M_fixed,"gray","InputType","bit");

        %small scale fading
        if isLOS
            %rician fading for LOS + scattering added
            s = sqrt(K_LOS/(K_LOS+1));
            sigma_sc = sqrt(1/(2*(K_LOS+1)));
            h = s + (sigma_sc*(randn + 1j*randn));

        else
            %Rayleigh for NLOS 
            h = (1/sqrt(2))*(randn + 1j*randn);
        end
        
        % we apply receive scaling so the avg rx power becomes Pr_W
        rx_clean = sqrt(Pr_w)*h*sym_fixed; 

        % AWGN
        noise = sqrt(N0_W/2)*(randn(size(rx_clean))+ 1j*randn(size(rx_clean)));
        rx_noisy = rx_clean + noise;

        %channel equalisation
        eq = rx_noisy/(sqrt(Pr_w)*h);

        demod_fixed = qamdemod(eq,M_fixed,'OutputType','bit');
        %counting bit errors
        biterrors_fixed = biterrors_fixed + sum(demod_fixed ~= bits_fixed);
        totalbits_fixed = totalbits_fixed + length(bits_fixed);

        % AMC scheme employment
        bits_amc = randi([0 1], Nsym*k_amc,1);
        sym_amc = qammod(bits_amc,M_amc,'gray','InputType','bit');
        % using the same channel 

        if isLOS
            %rician fading for LOS + scattering added
            s = sqrt(K_LOS/(K_LOS+1));
            sigma_sc = sqrt(1/(2*(K_LOS+1)));
            h = s + (sigma_sc*(randn + 1j*randn));
        else
            %Rayleigh for NLOS 
            h = (1/sqrt(2))*(randn + 1j*randn);
        end
        rx_clean_amc = sqrt(Pr_w)*h*sym_amc;
        noise_amc = sqrt(N0_W/2)*(randn(size(rx_clean_amc))+1j*randn(size(rx_clean_amc)));
        rx_noisy_amc = rx_clean_amc + noise_amc;
        eq_amc = rx_noisy_amc/(sqrt(Pr_w)*h);
        demod_amc = qamdemod(eq_amc,M_amc,'OutputType','bit');
        biterrors_amc = biterrors_amc + sum(demod_amc ~= bits_amc);
        totalbits_amc = totalbits_amc + length(bits_amc);
 end
    avgRxSNR_db(ip) = snr_acc / nframes;
    ber_fixed(ip) = biterrors_fixed/totalbits_fixed;
    ber_amc(ip) = biterrors_amc/totalbits_amc;
    %throughput 
    Tframe =  Nsysm / BW;
    throughput_fixed(ip) = (Nsysm*bps(fixedmodulation))/Tframe;
    throughput_amc(ip) = ((1-ber_amc(ip))*(Nsysm*mean([4,16],'all')))/Tframe;
    fprintf('Pt=%2d dBm: avgRxSNR=%.2f dB, BER_fixed=%.3e, BER_AMC=%.3e\n', ...
        pt_dbm, avgRxSNR_db(ip), ber_fixed(ip), ber_amc(ip));
end

figure;
subplot(1,2,1);
semilogy(avgRxSNR_db, ber_fixed, '-o', avgRxSNR_db, ber_amc, '-s','LineWidth',1.5);
grid on; 
xlabel('Average Rx SNR (dB)'); 
ylabel('BER'); 
legend('Fixed (QPSK)','AMC (QPSK/16QAM)','Location','southwest');
title('SISO BER vs Avg Rx SNR — UMi');
subplot(1,2,2);
plot(avgRxSNR_db, throughput_fixed/1e3, '-o', avgRxSNR_db, throughput_amc/1e3, '-s','LineWidth',1.5);
grid on; 
xlabel('Average Rx SNR (dB)'); 
ylabel('Throughput (kbits/s)'); 
legend('Fixed','AMC');
title('Illustrative throughput (uncoded)');


