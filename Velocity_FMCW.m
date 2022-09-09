close all;
clear all;
clc
tic
[y,fs] = audioread('COTS/testing/Range_Test_File.m4a');  % Load audio file (y: sampled data, fs: sample rate Hz)

% Seperate into sync data and radar backscatter data
data = y(:,1);          % radar backscatter
sync = y(:,2);          % sync data

f_start = 2.408e9;      % start frequency       Hz
f_stop = 2.495e9;       % stop frequency        Hz
BW = f_stop - f_start;  % Bandwidth Hz

Tp = 20e-3;             % Pulse width           s
N = Tp*fs;              % Number of samples

c = 2.997e8;            % speed of light        m/s

rr = c/(2*BW);          % range resolution      m
Rmax = N*rr/2;          % max unambiguous range m

%This is a change to test my git setup

