
clc; clear;

file = './Database/test_babble.wav';
[x, fs] = audioread(file);
len = length(x);
result = vad(x,fs,len);
figure(1);
subplot(2,1,1),plot(x),title('Signal'),axis([1 len min(x) max(x)]);
subplot(2,1,2),plot(result),title('VAD'),axis([1 len 0 1.2]);
