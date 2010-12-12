% CRANK
%
%   Demo over het krukdrijfstangmechanisme
%   Ter ondersteuning van het vak KINEMATICA EN WERKTUIGENDYNAMICA
% 
%   Copyright (c) 1999 by KULeuven, PMA

echo off
echo on, clc
% KRUK-DRIJFSTANGMECHANISME
% Dit programma illustreert de beweging van het kruk-drijfstangmechanisme
%
% Voor de demo kan je de waarden van de lengten van kruk en drijfstang zelf kiezen. 
% Voor het gemak is een standaardwaarde gegeven tussen haakjes. Druk op Enter of Return 
% om deze standaardwaarden te aanvaarden, of typ je eigen keuze.
%
echo off
fprintf('\nDruk op een toets om verder te gaan ...'), pause, disp(' ')
%
% Initialisatie van de parameters
R = input('Lengte R van de kruk (een positief getal vb 1)');
L = input('Lengte L van de drijfstang (L>R) vb 3');

N=1000

theta=(0:1:(N-1))*2*pi/N;

x=R*(1-cos(theta))+L*(1-sqrt(1-(R*sin(theta)/L).^2));
clf, hold off
plot(theta,x);
xlabel('\theta [rad]')
ylabel('Positie x')
title('Exacte oplossing voor de positie x')

disp('De eerste figuur toont de exacte oplossing voor de positie x voor 1 periode')
fprintf('\nDruk op een toets om verder te gaan ...\n'), pause, disp(' ')


X=fft(x,N);
a=2*real(X)/N;
b=-2*imag(X)/N;

clc
disp('We nemen de Fouriertransformatie van de exacte oplossing en bekijken')
disp('de harmonische inhoud aan de hand van de coeffienten van de Fourierreeks.')
format long
disp('De constante term')
a0 = a(1)/2
fprintf('\nDruk op een toets om verder te gaan ...\n'), pause, disp(' ')

disp('De coefficienten van de cosinus in de Fourierreeks:')
ak = a(2:11)'
disp('Merk op dat de coefficienten horende bij de hogere oneven harmonischen nul zijn.')
fprintf('\nDruk op een toets om verder te gaan ...\n'), pause, disp(' ')

disp('De coefficienten van de sinus in de Fourierreeks:')
bk = b(2:11)'
disp('Deze zijn niet exact nul, maar kunnen als nul beschouwd worden gezien de machinenauwkeurigheid.')
fprintf('\nDruk op een toets om verder te gaan ...\n'), pause, disp(' ')

X=[X(1) 2*X(2:floor(N/2))]/N;
figure
stem(0:10,abs(X(1:11)), 'o')
xlabel('Harmonische')
ylabel('Amplitude')
title('Frequentiespectrum')
zoom on
disp('De figuur toont het frequentiespectrum waarop te zien is dat de hogere')
disp('oneven harmonische nul zijn.')
fprintf('\nDruk op een toets om verder te gaan ...'), pause, disp(' ')

%
x2=a(1)/2+a(2)*cos(theta);
figure
plot(theta,x-x2)
xlabel('\theta [rad]')
ylabel('Amplitude')
title('Na verwijderen van DC-component en eerste harmonische')

%clc
disp('We nemen de constante term en de eerste harmonische van de reeksontwikkeling')
disp('en trekken deze af van de exacte vergelijking. Nu wordt de tweede harmonische')
disp('overheersend, maar ze heeft duidelijk een kleinere amplitude.')
fprintf('\nDruk op een toets om verder te gaan ...'), pause, disp(' ')

%
x3=a(1)/2 +a(2)*cos(theta)+a(3)*cos(2*theta);
figure
plot(theta,x-x3)
xlabel('\theta [rad]')
ylabel('Amplitude')
title('Na verwijderen van de tweede harmonische')

%clc
disp('We trekken nu ook de tweede harmonische af van dit resultaat waarna de vierde')
disp('harmonische (en niet de derde!) dominant wordt met nog kleinere amplitude.')
fprintf('\nDruk op een toets om verder te gaan ...\n'), pause, disp(' ')

%
x4=a(1)/2 +a(2)*cos(theta)+a(3)*cos(2*theta)+a(5)*cos(4*theta);
figure
plot(theta,x-x4)
xlabel('\theta [rad]')
ylabel('Amplitude')
title('Na verwijderen van de vierde harmonische')

%clc
disp('Analoog voor het verwijderen van de vierde harmonische: de zesde harmonische wordt zichtbaar.')
fprintf('\nDruk op een toets om verder te gaan ...\n'), pause, disp(' ')

%
x5=a(1)/2 +a(2)*cos(theta)+a(3)*cos(2*theta)+a(5)*cos(4*theta)+a(7)*cos(6*theta);
figure
plot(theta,x-x5)
xlabel('\theta [rad]')
ylabel('Amplitude')
title('Na verwijderen van de zesde harmonische')

%clc
disp('Tenslotte verwijderen we de zesde harmonische en de achtste harmonische wordt zichtbaar.')
