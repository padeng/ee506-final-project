% Mohammad Ismail Hossain 
% Jacobs University Bremen

% By using BPSK Constellation, Matlab simulation for SER of 1-by-1 system,
% Receive- MRC for 1-by-2 system, Transmit Beamforming and Alamouti Code 
% for 2-by-1 system. Channels and noises are all zero means with unit
% variances. 

clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 10^6; 
r_n = rand(1,N)>0.5;
BPSK = 2*r_n-1;
E_n_dB = -2:1:15; 
E_n_ln = 10.^(E_n_dB/10);
n_rcv =  [1 2];
n_trans=2;
th_Ber_ray = 0.5.*(1-sqrt(E_n_ln./(E_n_ln+1))); 
p_R_MRC = 1/2 - 1/2*(1+1./E_n_ln).^(-1/2);
th_Ber_R_MRC = p_R_MRC.^2.*(1+2*(1-p_R_MRC)); 
p_Alamouti = 1/2 - 1/2*(1+2./E_n_ln).^(-1/2);
th_Ber_Alamouti = p_Alamouti.^2.*(1+2*(1-p_Alamouti)); 
%%%%%%%%%%%%%%%%%%%%%%%%% One by One System  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:length(E_n_dB)
   W_n = 1/sqrt(2)*[randn(1,N) + 1j*randn(1,N)]; 
   C_h = 1/sqrt(2)*[randn(1,N) + 1j*randn(1,N)]; 
   R_x= C_h.*BPSK + 10^(-E_n_dB(p)/20)*W_n; 
   E_x = R_x./C_h;
   D_x = real(E_x)>0;
   N_Errors_ray(p) = size(find([r_n- D_x]),2);
end
S_Ber_ray = N_Errors_ray/N;
semilogy(E_n_dB,th_Ber_ray,'-g','LineWidth',2);hold on;
semilogy(E_n_dB,S_Ber_ray,'-pr','LineWidth',2);
%%%%%%%%%%%%%%%%%%% Receive MRC one by Two System  %%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:length(n_rcv)
    for q = 1:length(E_n_dB)
        W_n = 1/sqrt(2)*[randn(n_rcv(p),N) + 1j*randn(n_rcv(p),N)];
        C_h = 1/sqrt(2)*[randn(n_rcv(p),N) + 1j*randn(n_rcv(p),N)]; 
        M_BPSK = kron(ones(n_rcv(p),1),BPSK);
        R_x = C_h.*M_BPSK + 10^(-E_n_dB(q)/20)*W_n;
        E_x = sum(conj(C_h).*R_x,1)./sum(C_h.*conj(C_h),1); 
        D_x = real(E_x)>0;
        N_Errors_MRC(p,q) = size(find([r_n- D_x]),2);
    end
end
S_Ber_R = N_Errors_MRC/N;
S_Ber_R_MRC = S_Ber_R(2,:);
semilogy(E_n_dB,th_Ber_R_MRC,'-r','LineWidth',2);hold on;
semilogy(E_n_dB,S_Ber_R_MRC,'-ob','LineWidth',2);
%%%%%%%%%%%%%% Transmit beamforming for two by One System  %%%%%%%%%%%%%%%%
for p = 1:length(E_n_dB)
    W_n = 1/sqrt(2)*[randn(1,N) + 1j*randn(1,N)]; 
    C_h = 1/sqrt(2)*[randn(n_trans,N) + 1j*randn(n_trans,N)]; 
    M_BPSK = (1/sqrt(n_trans))*kron(ones(n_trans,1),BPSK); 
    C_h_E = C_h.*exp(-1j*angle(C_h));
    R_x = sum(C_h_E.*M_BPSK,1) + 10^(-E_n_dB(p)/20)*W_n; 
    E_x = R_x./sum(C_h_E,1); 
    D_x = real(E_x)>0;
    N_Errors_Beam(p) = size(find([r_n- D_x]),2);
end
S_Ber_T_Beam = N_Errors_Beam/N; 
semilogy(E_n_dB,S_Ber_T_Beam,'-*m','LineWidth',2);
%%%%%%%%%%%%%%%%% Alamouti Code for two by One System  %%%%%%%%%%%%%%%%%%%%
for p = 1:length(E_n_dB)
    Z_c = zeros(2,N);
    Z_c(:,1:2:end) = (1/sqrt(2))*reshape(BPSK,2,N/2);
    P_c=flipud(reshape(conj(BPSK),2,N/2));
    Z_c(:,2:2:end) = (1/sqrt(2))*(kron(ones(1,N/2),[-1;1]).*P_c); 
    C_h = 1/sqrt(2)*[randn(1,N) + 1j*randn(1,N)];
    M_BPSK = kron(reshape(C_h,2,N/2),ones(1,2));     
    W_n = 1/sqrt(2)*[randn(1,N) + 1j*randn(1,N)];
    R_x = sum(M_BPSK.*Z_c,1) + 10^(-E_n_dB(p)/20)*W_n;
    M_Rx = kron(reshape(R_x,2,N/2),ones(1,2)); 
    M_Rx(2,:) = conj(M_Rx(2,:)); 
    Z_E = zeros(2,N);
    Z_E(:,[1:2:end]) = reshape(C_h,2,N/2); 
    Z_c=flipud(reshape(C_h,2,N/2));
    Z_E(:,[2:2:end]) = kron(ones(1,N/2),[1;-1]).*Z_c; 
    Z_E(1,:) = conj(Z_E(1,:));
    Z_E_P = sum(Z_E.*conj(Z_E),1);
    E_x = sum(Z_E.*M_Rx,1)./Z_E_P; 
    E_x(2:2:end) = conj(E_x(2:2:end));
    D_x = real(E_x)>0;
    N_Errors_Alamouti(p) = size(find([r_n- D_x]),2);
end
S_Ber_Alamouti = N_Errors_Alamouti/N; 
semilogy(E_n_dB,th_Ber_Alamouti,'-+c','LineWidth',2);hold on;
semilogy(E_n_dB,S_Ber_Alamouti,'-dk','LineWidth',2);
grid on
legend('Theory-(1-by-1) System', 'Simulation-(1-by-1) System','Theory-Rx-MRC (1-by-2) System', 'Simulation-Rx-MRC (1-by-2) System','Simulation-Tx-Beam (2-by-1) System','Theory-Alamouti (2-by-1) System', 'Simulation-Alamouti (2-by-1) System','');
xlabel('SNR (dB)','fontsize',14);
ylabel('Symbol Error Rate (SER)','fontsize',14);
title('SER for BPSK modulation','fontsize',14);
axis([-2 max(E_n_dB) min(S_Ber_R_MRC) max(S_Ber_ray)])