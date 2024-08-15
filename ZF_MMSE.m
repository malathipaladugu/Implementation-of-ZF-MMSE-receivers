close all;
clear all;

block_len = 1000;
n_blocks = 10000;
t=2 ; r=2;% #rx_ant = tx_ant =2


Eb_dB = 1.0:4.0:33.0;%signal power varied from 1dB to 33dB
Eb = 10.^(Eb_dB/10);%Signal power in dB
No = 1;%Noise power = 1
Es = 2*Eb;%Symbol energy
SNR = Es/No;
SNR_dB = 10*log10(SNR);

%To caluclate BER at different signal powers
BER_ZF = zeros(1,length(Eb_dB)); %For ZF receiver
BER_LMMSE = zeros(1,length(Eb_dB));

%ZF_receiver x_hat = pinv(H)*Y
%MIMO System model is Y= HX +V
for blk = 1:n_blocks
    %Defining tansmitted symbol X and converting chunks to symbol(QPSK)
    I_bits = randi([0,1],t,block_len);
    Q_bits = randi([0,1],t,block_len);
    X = (2*I_bits-1)+1j*(2*Q_bits-1);
    %Defining channel coefficient matrix H
    H = 1/sqrt(2)*(randn(r,t)+ 1j*randn(r,t));
    %Defining the noise samples V
    V = sqrt(No/2)*(randn(r, block_len)+ 1j* randn(r, block_len));

    for k = 1:length(Eb_dB)
        X = sqrt(Eb(k))*X;
        Y = H*X + V;

        %Calculating the ZF receiver estimate X_hat_ZF = pinv(H)*Y
        X_hat_ZF = pinv(H)*Y;
        %Decoding the received bits
        %(1,1)= (1+j), (0,0) = (-1-j), (0,1)= (-1+j), (1,0) = (1-j)
        Dec_I_bits_ZF = (real(X_hat_ZF)>0); 
        Dec_Q_bits_ZF = (imag(X_hat_ZF)>0);
        %BER calculation by comparing the transmitted and decoded symbols
        BER_ZF(k) = BER_ZF(k) + sum(sum(Dec_I_bits_ZF ~= I_bits)) + sum(sum(Dec_Q_bits_ZF ~= Q_bits));

        %calculating LMMSE receiver estimate X_hat_LMMSE =
        %inv((transpose(H)*H + I/SNR)*transpose(H)*Y
        X_hat_LMMSE = inv(H'*H+ eye(t)/(Es(k)/No))* H' * Y;
        %Decoding the received bits
        Dec_I_bits_LMMSE = (real(X_hat_LMMSE)>0); 
        Dec_Q_bits_LMMSE = (imag(X_hat_LMMSE)>0);
        %BER calculation for LMMSE estimate
        BER_LMMSE(k) = BER_LMMSE(k) + sum(sum(Dec_I_bits_LMMSE ~= I_bits)) + sum(sum(Dec_Q_bits_LMMSE ~= Q_bits));
    end
end

BER_ZF = BER_ZF/block_len/n_blocks/2/t;
BER_LMMSE = BER_LMMSE/block_len/n_blocks/2/t;
L = r-t+1;
BERt = nchoosek(2*L-1, L)/2^L./SNR.^L;

semilogy(SNR_dB, BER_ZF, 'b s', 'linewidth',3.0, 'MarkerFaceColor','b','MarkerSize',9.0);
hold on;
semilogy(SNR_dB, BER_LMMSE,'r -- o','linewidth',3.0, 'MarkerFaceColor','r','MarkerSize',9.0);
semilogy(SNR_dB, BERt, 'g -.', 'linewidth',3.0, 'MarkerFaceColor','g','MarkerSize',9.0);
axis tight;
grid on;
legend('ZF','LMMSE', 'ZF_Theoretical');
xlabel('SNR(dB)');
ylabel('BER');
title('BER vs SNR(dB) for MIMO ZF, LMMSE receivers');


    