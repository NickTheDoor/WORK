function [] = nmesfunction()

%This routine works for Matheus project
%by fernando - fdiefenthaeler@gmail.com
%09/10/2017
%by nicola - nicolalaporta9@gmail.com
%~~/11/2017
%=========================================================================
clear
close all
clc
%% import .edf RTD file
[filename,pathname] = uigetfile('*.edf');
%          if exist([filename(1:end-4),'.mat'])
%          load([filename(1:end-4),'.mat'])
%          else
fulladdress = [pathname filename];

Data = importdata(fulladdress);
[hdr, data] = edfRead((fulladdress));

%% name each data
gain1 = 100;
ta = (detrend(data(1,:))/gain1)';
sinc = data(2,:)';
torque = data(3,:)';

%% torque 
figure(1);
plot(torque)
fator = -1; %input('digite -1 para girar o gráfico 180º = ');

close

%% filtra torque
SR = 2000;        %sample frequency
SRr = 2*(10/SR);
[B,A] = butter(5,SRr,'low');
Torquefilt = filtfilt(B,A,torque)*fator;
    %Pega o menor valor do vetor Torquefilt
    menorVal = abs(min(Torquefilt));
clear('A','B')
%figure(2);
   
    figure('units','normalized','outerposition',[0 0 1 1])
plot(ta,'r')
    hold on
plot(Torquefilt,'k')

question = {'O gráfico apresenta CVIM e twitches pós CVIM? [s/n]'};
titulo = 'Observe o gráfico';
dims = [1 50];
definput = {'s'};
opts.Resize = 'on';
answer = inputdlg(question,titulo,dims,definput,opts);
resposta = answer{:};

if resposta == 's' 
    uiwait(msgbox({'Os dois(2) primeiros serão para a baseline do sinal.' ''...
 'Sendo que o segundo (2º) deverá ser no começo das 250 contrações.' 'O terceiro (3º) ao final das 250 contrações.'...
 'E o quarto (4º) anterior ao começo da CVIM.'},'Selecione os pontos','modal'));

    uiwait(msgbox({ 'ESPAÇO para começar a seleção, ENTER para finalizar' ''...
 'Em caso de seleção errada, BACKSPACE para deletar o ultimo ponto selecionado'},'Selecione os 4 pontos','modal'));


if exist([filename(1:end-4) ' ponto1.mat']) == 2
    
    load ([filename(1:end-4) ' ponto1.mat'])
    
     uiwait(msgbox({ 'O ponto ja foi selecionado previamente'},'modal'));

    
else
[pontos, valory] =  ginput_zoom(4);
    ponto = round(pontos);
    close all
    save([filename(1:end-4) ' ponto1.mat'],'ponto')      
end
    % todas as funções de calculo
    
elseif resposta == 'n'
    uiwait(msgbox({'Os dois(2) primeiros serão para a baseline do sinal.' ''...
 'Sendo que o segundo (2º) deverá ser no começo das 250 contrações.' '' 'O terceiro (3º) ao final das 250 contrações.'...
 '' 'E o quarto (4º) Após os ultimos twitches'},...
 'Selecione os pontos','modal'));

    uiwait(msgbox({ 'ESPAÇO para começar a seleção, ENTER para selecionar' ''...
 'Em caso de seleção errada, BACKSPACE para deletar o ultimo ponto selecionado'},'Selecione os 4 pontos','modal'));

if exist([filename(1:end-4) ' ponto1.mat']) == 2
    
    load ([filename(1:end-4) ' ponto1.mat'])
    
      uiwait(msgbox({ 'O ponto ja foi selecionado previamente'},'modal'));
      
else
[pontos, valory] =  ginput_zoom(4);
    ponto = round(pontos);
    close all
    save([filename(1:end-4) ' ponto1.mat'],'ponto')      
end
    
    %fazer funções de calculo até antes da CVIM
    
else    
    uiwait(msgbox({'Escreva uma opção válida'},'ERRO','modal'));
end


% [pontos, valory] =  ginput_zoom(4);
%     ponto = round(pontos);

 torqueCortado = Torquefilt(ponto(2):ponto(3));

%% Torque cortado positivo
 torqueCortadoPos = torqueCortado + menorVal;
  
%% Tara EMG
mediaEMG = mean(ta(ponto(1):ponto(2)));

taTara = ta + abs(mediaEMG);
 
%% Tara Torquefilt
 mediaTq = mean(Torquefilt(ponto(1):ponto(2))); 
 TqTara = Torquefilt + abs(mediaTq); 

close all;

%% Plota e calcula a média do plato CVM

if resposta == 's' 

if exist([filename(1:end-4) ' platoCVIM.mat']) == 2
    load ([filename(1:end-4) ' platoCVIM.mat'])
     uiwait(msgbox({'O ponto ja foi selecionado previamente'},'Selecione 1 ponto','modal'));
    close all
else
   
    figure('units','normalized','outerposition',[0 0 1 1])  
plot(Torquefilt(ponto(4):ponto(4)+15000),'k','LineWidth',1);
title('CVIM');
    uiwait(msgbox({'Selecione um ponto para determinar a média do platô de torque.' ''...
'ESPAÇO para começar a seleção, ENTER para selecionar' ''...
'Em caso de seleção errada, BACKSPACE para deletar o ultimo ponto selecionado'},'Selecione 1 ponto','modal'));
   
[platox, ~] = ginput_zoom(1);
platox1 = round(platox + ponto(4));   
    close all
    save([filename(1:end-4) ' platoCVIM.mat'],'platox1')      
end
    
   
    mediaPlato = mean(TqTara(platox1-400:platox1+500));
   
%     mediaPlato = mean(TqTara(platox-500:platox+500));
    
    
% platox1 = round(platox(1) + ponto(4));
% platox2 = round(platox(2) + ponto(4));
%mediaPlato = mean(Torquefilt(platox1:platox2));
%mediaPlato = mean(Torquefilt(platox(1)+ ponto(4):platox(2)+ ponto(4)));

close all;

end

% %% Calcula média da m-wave durante platô da CVIM
% 
% plot(TqTara(ponto(4):ponto(4)+22500),'k','LineWidth',1);
% hold on
% plot(taTara(ponto(4):ponto(4)+22500),'k','LineWidth',1);
% hold off
% title('CVIM');
% 
% periodoMwave = (platox1-500:platox1+500);
% valoresMwave = taTara(platox1-500:platox1+500);
% mediaMwaveCVIM = rms(valoresMwave);
% 
% fprintf('A média da eletromiografia durante o platô da CVIM é: %.2f mV \n', mediaMwaveCVIM);


%% Algoritmo que seleciona m-waves pré protocolo e as analiza

%       function [duracaoTq1 amplitudeTq1 timetoPeakTq1  halfRelaxTimeTq1  tdf1]...
%                = twitches1(taTara, TqTara, ponto);
  

                                               
    [mPicos mLocais] = findpeaks(taTara(1:ponto(2)),'MinPeakHeight',4,'MinPeakDistance',1000);   
    [tPicos tLocais] = findpeaks(TqTara(1:ponto(2)),'MinPeakProminence',0.95,'MinPeakDistance',1000);      
    
%Onsets dos twitchs
    [mPicos1 mLocais1] = findpeaks(taTara(1:ponto(2))*-1,'MinPeakHeight',4,'MinPeakDistance',1000);   
onsetM1 = mPicos1 * -1;    

% Acha o onset plotando o gráfico de torque dentro do intervalo 1:ponto(2) 
% adicionando +10 ao torque para deixa-lo positivo

graficoT1 = TqTara(1:ponto(2))+10;  
localOnsetT1 = tLocais;
     
for i = 1:3
     while abs(graficoT1(localOnsetT1(i) - 15)) - abs(graficoT1 (localOnsetT1(i) - 25)) > 0.001
        localOnsetT1(i) = localOnsetT1(i) - 10;     
     end
end

onsetT1 = graficoT1(localOnsetT1(1:3))-10;

% Achar tempo de meio relaxamento
  
  graficoMeio = graficoT1;
    
amplitudeT1 = tPicos(1:3) - onsetT1(1:3);
    valorInicio1 = tPicos(1:3);
valorMeio1 = (tPicos(1:3)-(amplitudeT1/2))+10;   

    valorFim = onsetT1;
localMeio1 = tLocais(1:3);
    
for i1 = 1:3
    while graficoMeio(localMeio1(i1)) > valorMeio1(i1)   
        localMeio1(i1) = localMeio1(i1) + 1;
    end
end    
    valorMeio1 = valorMeio1 - 10;       
      
 %Achar o local de termino da contração 
    localContinua1 = localMeio1;        
for i11 = 1:3
%  zero = find(~pontoCritico); 
derivada1 = diff(graficoMeio(localContinua1(i11):localContinua1(i11)+500));
 
    [pontoCritico1(i11), localCrit1(i11)] = min(abs(derivada1));    
    
localCrit1(i11) = localCrit1(i11) + localContinua1(i11);
end   
    pontoCritico1 = graficoMeio(localCrit1) - 10;
 
% Plot das m-waves pré protocolo (eletroestimulação + torque)
    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    addToolbarExplorationButtons(fig1);
    plot(taTara(1:ponto(2)), 'r','LineWidth',0.2);
hold on   
    plot(TqTara(1:ponto(2)),'k','LineWidth',1);
    plot(mLocais(1:3),mPicos(1:3),'bo');
    plot(tLocais(1:3),tPicos(1:3),'b^'); 
    plot(localOnsetT1(1:3), onsetT1(1:3),'bv');
    plot(mLocais1(1:3), onsetM1(1:3),'bo');
    plot(localMeio1, valorMeio1,'c*'); 
    plot(localCrit1,pontoCritico1,'bv');
hold off
    title('Sinal pré protocolo');
    
 
    
    question = {[' Se sim, digite: 0.' sprintf('\n')...
     ' Caso não estejam, qual deve ser a nova amplitude mínima para que esses picos sejam incluidos na função de busca?'...
     sprintf('\n') 'Amplitude atual: 0.5 ' sprintf('\n') 'Ex: 0.7       Ex2: 0.4']};
    titulo = 'Os picos de torque se encontram no lugar certo?';
    dims = [1 70];
    definput = {'0'};
%   opts.Resize = 'on';
    answer1 = inputdlg(question,titulo,dims,definput); %opts);
    resposta1 = str2num(answer1{:});

    %Função para checar se os picos de torque foram selecionados
    %corretamente

  while resposta1 ~= 0
    
 close;
 
 % Plot para analisar as m-waves pré protocolo (eletroestimulação + torque)
    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    addToolbarExplorationButtons(fig1);
   
hold on   
    plot(TqTara(1:ponto(2)),'k','LineWidth',1);
    plot(mLocais(1:3),mPicos(1:3),'bo');
    plot(tLocais(1:3),tPicos(1:3),'b^'); 
    plot(localOnsetT1(1:3), onsetT1(1:3),'bv');
    plot(mLocais1(1:3), onsetM1(1:3),'bo');
    plot(localMeio1, valorMeio1,'c*'); 
    plot(localCrit1,pontoCritico1,'bv');
hold off
    title('Sinal pré protocolo'); 
 
        question = {[' Se sim, digite: 0.' sprintf('\n')...
     ' Caso não estejam, qual deve ser a nova amplitude mínima para que esses picos sejam incluidos na função de busca?'...
     sprintf('\n') 'Amplitude atual: 0.5 ' sprintf('\n') 'Ex: 0.7       Ex2: 0.4']};
    titulo = ' E agora, os picos de torque se encontram no lugar certo?';
    dims = [1 70];
    definput = {'0'};
%   opts.Resize = 'on';
    answer1 = inputdlg(question,titulo,dims,definput) %opts);
    resposta1 = str2num(answer1{:});

%Repete a função de achar os picos de torque relacionados aos twitches      
    
    [mPicos mLocais] = findpeaks(taTara(1:ponto(2)),'MinPeakHeight',7,'MinPeakDistance',1000);   
    [tPicos tLocais] = findpeaks(TqTara(1:ponto(2)),'MinPeakProminence',resposta1,'MinPeakDistance',1000);      



%Onsets dos twitchs
    [mPicos1 mLocais1] = findpeaks(taTara(1:ponto(2))*-1,'MinPeakHeight',7,'MinPeakDistance',1000);   
onsetM1 = mPicos1 * -1;    

% Acha o onset de torque
graficoT1 = TqTara(1:ponto(2))+10;  
localOnsetT1 = tLocais;
     
for i = 1:3
     while abs(graficoT1(localOnsetT1(i) - 10)) - abs(graficoT1 (localOnsetT1(i) - 20)) > 0.001
        localOnsetT1(i) = localOnsetT1(i) - 10;     
     end
end


onsetT1 = graficoT1(localOnsetT1(1:3))-10;

% Achar tempo de meio relaxamento
  
  graficoMeio = graficoT1;
    
amplitudeT1 = tPicos(1:3) - onsetT1(1:3);
    valorInicio1 = tPicos(1:3);
valorMeio1 = (tPicos(1:3)-(amplitudeT1/2))+10;   

    valorFim = onsetT1;
localMeio1 = tLocais(1:3);
    
for i1 = 1:3
    while graficoMeio(localMeio1(i1)) > valorMeio1(i1)   
        localMeio1(i1) = localMeio1(i1) + 1;
    end
end    
    valorMeio1 = valorMeio1 - 10;       
      
 %Achar o local de termino da contração 
    localContinua1 = localMeio1;        
for i11 = 1:3
%  zero = find(~pontoCritico); 
derivada1 = diff(graficoMeio(localContinua1(i11):localContinua1(i11)+500));
 
    [pontoCritico1(i11), localCrit1(i11)] = min(abs(derivada1));    
    
localCrit1(i11) = localCrit1(i11) + localContinua1(i11);
end   
    pontoCritico1 = graficoMeio(localCrit1) - 10;
    
    end 
   
    %salva 1 twitches pré protocolo
    saveas(fig1, '1 twitches pré protocolo.fig')
    
    uiwait(msgbox('Será feita a análise dos primeiros twitchs','Twitchs','modal'));  
    
%Plot das 3 m-waves (uma embaixo da outra)    
    fig11 = figure('units','normalized','outerposition',[0 0 1 1]);
    addToolbarExplorationButtons(fig11);
subplot(3,1,1)
    plot(taTara(mLocais(1)-50:mLocais(1)+50), 'r','LineWidth',0.2);
    title('Análise individual das M-waves pré-protocolo, 1ª');  
subplot(3,1,2) 
    plot(taTara(mLocais(2)-50:mLocais(2)+50), 'r','LineWidth',0.2);
    title('Análise individual das M-waves pré-protocolo, 2ª');
subplot(3,1,3) 
    plot(taTara(mLocais(3)-50:mLocais(3)+50), 'r','LineWidth',0.2); 
    title('Análise individual das M-waves pré-protocolo, 3ª');    

   %salva 1 twitches pré protocolo, os 3
    saveas(fig11, '1 twitches pré protocolo, os 3.fig')  

%plot da primeira m-wave pré-protocolo  
fig12 = figure('units','normalized','outerposition',[0 0 1 1]);
addToolbarExplorationButtons(fig12);
plot(taTara(mLocais(1)-50:mLocais(1)+90), 'r','LineWidth',0.2);
title('Análise individual das M-waves pré-protocolo, 1ª');
    
   %salva 1 m-waves pré protocolo - primeira
   saveas(fig12, '1 m-waves pré protocolo - primeira.fig')  

%plot da segunda m-wave pré-protocolo  
fig13 = figure('units','normalized','outerposition',[0 0 1 1]);
addToolbarExplorationButtons(fig13);
plot(taTara(mLocais(2)-50:mLocais(2)+90), 'r','LineWidth',0.2); 
title('Análise individual das M-waves pré-protocolo, 2ª'); 

 %salva 1 m-waves pré protocolo - segunda
   saveas(fig13, '1 m-waves pré protocolo - segunda.fig')  

%plot da terceira m-wave pré-protocolo     
fig14 = figure('units','normalized','outerposition',[0 0 1 1]);
addToolbarExplorationButtons(fig14);
plot(taTara(mLocais(3)-50:mLocais(3)+90), 'r','LineWidth',0.2); 
title('Análise individual das M-waves pré-protocolo, 3ª');    
  
 %salva 1 m-waves pré protocolo - terceira
   saveas(fig14, '1 m-waves pré protocolo - terceira.fig')  

%  uiwait(msgbox({'A análise das m-waves é avaliador-dependente' ''...
%     'Utilize o "Data cursor" para a análise' '' 'Botão com uma caixinha amarela e uma curva passando por ela'}...
%     ,'Análise m-waves','modal'));

%Cálculo e armazenamento das variáveis relevantes
 duracaoTq1 = localCrit1' - localOnsetT1(1:3);
 amplitudeTq1 = tPicos(1:3) - onsetT1(1:3);
 timetoPeakTq1 = tLocais(1:3) - localOnsetT1(1:3);
 halfRelaxTimeTq1 = localMeio1 - tLocais(1:3);
 tdf1 = amplitudeTq1 ./ timetoPeakTq1;   

%     end
 %close all
%% Algoritmo que seleciona m-waves pré CVIM e as analiza
%essa é a tentativa q deu realmente certo
                                                                       %'MinPeakHeight',7,%             
[mPicos2 mLocais2] = findpeaks(taTara(ponto(3):ponto(4)),'MinPeakHeight',5,'MinPeakDistance',2000);  

[tPicos2 tLocais2] = findpeaks(TqTara(ponto(3):ponto(4)),'MinPeakProminence',0.5,'MinPeakDistance',1000);   

graficoMwave2 = taTara(ponto(3):ponto(4));
graficoTorque2 = TqTara(ponto(3):ponto(4));


% xx1 = 2952;
% xx2 = 9942;
% xx3 = 18830
% 
% tPicos2(1) = graficoTorque2(xx1);
% tLocais2(1) = xx1;
% 
% tPicos2(2) = graficoTorque2(xx2);
% tLocais2(2) = xx2;
% 
% tPicos2(3) = graficoTorque2(xx3);
% tLocais2(3) = xx3;

% tGraf = TqTara(ponto(3):ponto(4));
% gPeaks = tGraf(mLocais2 -1500 : mLocais2 +1500);
% 
% [tPicos2 tLocais2] = findpeaks(gPeaks),'MinPeakProminence',0.5,'MinPeakDistance',500;   




%Onsets dos twitchs
[mPicos22 mLocais22] = findpeaks(taTara(ponto(3):ponto(4))*-1,'MinPeakHeight',7,'MinPeakDistance',2000);   
    onsetM2 = mPicos22 *-1;
    
%  Acha o onset plotando o gráfico de torque dentro do intervalo ponto(3) ao
% ponto(4) adicionando +10 ao torque para deixa-lo positivo e achar um ponto
% crítico na curva e ao final subtrai esses +10 para normalizar


graficoT2 = TqTara(ponto(3):ponto(4))+10;
localOnsetT2 = tLocais2;

% localOnsetT2 = [tLocais2(1), tLocais2(3), tLocais2(5)]; %tLocais2;
     
for ii = 1:3
     while abs(graficoT2(localOnsetT2(ii)-15)) - abs(graficoT2 (localOnsetT2(ii) - 25)) > 0.001
localOnsetT2(ii) = localOnsetT2(ii) - 10;     
     end
end
onsetT2 = graficoT2(localOnsetT2(1:3))-10;



%% Quando der erro

% fig999 = figure('units','normalized','outerposition',[0 0 1 1]);
% 
% plot(TqTara(ponto(3):ponto(4)),'k','LineWidth',1);
% hold on
% plot(taTara(ponto(3):ponto(4)),'r','LineWidth',1);
% hold off
%%
% Achar tempo de meio relaxamento
%     graficoMeio2 = graficoT2-10;
    
    graficoMeio2 = graficoT2;
    
amplitudeT2 = tPicos2(1:3) - onsetT2(1:3)';
    valorInicio2 = tPicos2(1:3);
valorMeio2 = (tPicos2(1:3)-(amplitudeT2/2))+10;   

    valorFim2 = onsetT2;
localMeio2 = tLocais2(1:3);%[tLocais2(1), tLocais2(3), tLocais2(5)]; %tLocais2(1:3);
    
for i2 = 1:3
    while graficoMeio2(localMeio2(i2)) > valorMeio2(i2)   
        localMeio2(i2) = localMeio2(i2) + 1;
    end
end
    
valorMeio2 = valorMeio2 - 10;   
       

%Achar o local de termino da contração 
localContinua2 = localMeio2;        

for i22 = 1:3
%  zero = find(~pontoCritico); 
derivada2 = diff(graficoMeio2(localContinua2(i22):localContinua2(i22)+500));
 
    [pontoCritico2(i22), localCrit2(i22)] = min(abs(derivada2));    
    
localCrit2(i22) = localCrit2(i22) + localContinua2(i22);

end   
pontoCritico2 = graficoMeio2(localCrit2) - 10;

 % Plot para analisar as m-waves pré CVIM (eletromiografia + torque)
fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
addToolbarExplorationButtons(fig2);
    if resposta == 's' 
 plot(taTara(ponto(3):ponto(4)+10000), 'r','LineWidth',0.2);
    hold on   
 plot(TqTara(ponto(3):ponto(4)+10000),'k','LineWidth',1);
    elseif resposta == 'n'
 plot(taTara(ponto(3):ponto(4)), 'r','LineWidth',0.2);
    hold on   
 plot(TqTara(ponto(3):ponto(4)),'k','LineWidth',1);        
    end  
 plot(mLocais2(1:3), mPicos2(1:3),'bo');  
plot(tLocais2(1:3), tPicos2(1:3),'b^');  %  plot([tLocais2(1), tLocais2(3), tLocais2(5)],[tPicos2(1), tPicos2(3), tPicos2(5)], 'b^'); 
 plot(localOnsetT2(1:3), onsetT2(1:3) ,'bv');
 plot(mLocais22(1:3), onsetM2(1:3),'bo')
 plot(localMeio2, valorMeio2,'c*'); 
 plot(localCrit2,pontoCritico2,'bv');
    hold off 
    title('Sinal pré CVIM');
 
            
  question = {[' Se sim, digite: 0.' sprintf('\n')...
     ' Caso não estejam, qual deve ser a nova amplitude mínima para que esses picos sejam incluidos na função de busca?'...
     sprintf('\n') 'Amplitude atual: 0.5 ' sprintf('\n') 'Ex: 0.7       Ex2: 0.4']};
    titulo = 'Os picos de torque se encontram no lugar certo?';
    dims = [1 70];
    definput = {'0'};
    answer2 = inputdlg(question,titulo,dims,definput);
    resposta2 = str2num(answer2{:});

    %Função para checar se os picos de torque foram selecionados
    %corretamente

   while resposta2 ~= 0
    
 close;   
    
 % Plot das 3 m-waves pré CVIM
fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
addToolbarExplorationButtons(fig2);
    if resposta == 's' 
 plot(taTara(ponto(3):ponto(4)+10000), 'r','LineWidth',0.2);
    hold on   
 plot(TqTara(ponto(3):ponto(4)+10000),'k','LineWidth',1);
    elseif resposta == 'n'
 plot(taTara(ponto(3):ponto(4)), 'r','LineWidth',0.2);
    hold on   
 plot(TqTara(ponto(3):ponto(4)),'k','LineWidth',1);        
    end  
 plot(mLocais2(1:3), mPicos2(1:3),'bo');  
 plot(tLocais2(1:3), tPicos2(1:3),'b^');  
 plot(localOnsetT2(1:3), onsetT2(1:3) ,'bv');
 plot(mLocais22(1:3), onsetM2(1:3),'bo')
 plot(localMeio2, valorMeio2,'c*'); 
 plot(localCrit2,pontoCritico2,'bv');
    hold off 
    title('Sinal pré CVIM');
     
  question = {[' Se sim, digite: 0.' sprintf('\n')...
     ' Caso não estejam, qual deve ser a nova amplitude mínima para que esses picos sejam incluidos na função de busca?'...
     sprintf('\n') 'Amplitude atual: 0.5 ' sprintf('\n') 'Ex: 0.7       Ex2: 0.4']};
    titulo = ' E agora, os picos de torque se encontram no lugar certo?';
    dims = [1 70];
    definput = {'0'};
    answer2 = inputdlg(question,titulo,dims,definput);
    resposta2 = str2num(answer2{:});    
    
%Repete a função de achar os picos de torque relacionados aos twitches  
    
[mPicos2 mLocais2] = findpeaks(taTara(ponto(3):ponto(4)),'MinPeakHeight',7,'MinPeakDistance',2500);   
[tPicos2 tLocais2] = findpeaks(TqTara(ponto(3):ponto(4)),'MinPeakProminence',resposta2,'MinPeakDistance',2000);   

graficoMwave2 = taTara(ponto(3):ponto(4));


%Onsets dos twitchs
[mPicos22 mLocais22] = findpeaks(taTara(ponto(3):ponto(4))*-1,'MinPeakHeight',7,'MinPeakDistance',2000);   
    onsetM2 = mPicos22 *-1;
    
%  Acha o onset plotando o gráfico de torque dentro do intervalo ponto(3) ao
% ponto(4) adicionando +10 ao torque para deixa-lo positivo e achar um ponto
% crítico na curva e ao final subtrai esses +10 para normalizar

graficoT2 = TqTara(ponto(3):ponto(4))+10;
localOnsetT2 = tLocais2;
     
for ii = 1:3
     while abs(graficoT2(localOnsetT2(ii) - 15 ) - abs(graficoT2 (localOnsetT2(ii) - 25)) > 0.001)
localOnsetT2(ii) = localOnsetT2(ii) - 10;     
     end
end
onsetT2 = graficoT2(localOnsetT2(1:3))-10;

% Achar tempo de meio relaxamento
%     graficoMeio2 = graficoT2-10;
    
    graficoMeio2 = graficoT2;
    
amplitudeT2 = tPicos2(1:3) - onsetT2(1:3);
    valorInicio2 = tPicos2(1:3);
valorMeio2 = (tPicos2(1:3)-(amplitudeT2/2))+10;   

    valorFim2 = onsetT2;
localMeio2 = tLocais2(1:3);
    
for i2 = 1:3
    while graficoMeio2(localMeio2(i2)) > valorMeio2(i2)   
        localMeio2(i2) = localMeio2(i2) + 1;
    end
end
    
valorMeio2 = valorMeio2 - 10;   
       

%Achar o local de termino da contração 
localContinua2 = localMeio2;        

for i22 = 1:3
%  zero = find(~pontoCritico); 
derivada2 = diff(graficoMeio2(localContinua2(i22):localContinua2(i22)+500));
 
    [pontoCritico2(i22), localCrit2(i22)] = min(abs(derivada2));    
    
localCrit2(i22) = localCrit2(i22) + localContinua2(i22);

end   
pontoCritico2 = graficoMeio2(localCrit2) - 10; 
    
   end 
   
   %salva 2 twitches pós protocolo
   saveas(fig2, '2 twitches pós protocolo.fig')  
   
%     uiwait(msgbox('Será feita a análise dos Twitchs pré CVIM','Twitchs','modal'));



%plot das 3 m-waves
     fig21 = figure('units','normalized','outerposition',[0 0 1 1]);
     addToolbarExplorationButtons(fig21);
      
subplot(3,1,1)
    plot(graficoMwave2(mLocais2(1)-50:mLocais2(1)+90), 'r','LineWidth',0.2);
    title('Análise individual das M-waves pré-CVIM, 1ª');  
subplot(3,1,2)        
    plot(graficoMwave2(mLocais22(2)-50:mLocais22(2)+90), 'r','LineWidth',0.2);    
    title('Análise individual das M-waves pré-CVIM, 2ª');
subplot(3,1,3) 
    plot(graficoMwave2(mLocais2(3)-50:mLocais2(3)+90), 'r','LineWidth',0.2); 
    title('Análise individual das M-waves pré-CVIM, 3ª');    

    %salva 2 twitches pós protocolo, os 3
   saveas(fig21, '2 twitches pós protocolo, os 3.fig')  
close;    
% uiwait(msgbox({'A análise das m-waves é avaliador-dependente' ''...
%     'Utilize o "Data cursor" para a análise' '' 'Botão com uma caixinha amarela e uma curva passando por ela'}...
%     ,'Análise m-waves','modal'));  

fig22 = figure('units','normalized','outerposition',[0 0 1 1]);
addToolbarExplorationButtons(fig22);
 plot(graficoMwave2(mLocais2(1)-50:mLocais2(1)+90), 'r','LineWidth',0.2);
    title('Análise individual das M-waves pré-CVIM, 1ª');  
 
    %salva 2 m-waves pós protocolo - primeira
   saveas(fig22, '2 m-waves pós protocolo - primeira.fig')  
    
fig23 = figure('units','normalized','outerposition',[0 0 1 1]);
addToolbarExplorationButtons(fig23);
 plot(graficoMwave2(mLocais22(2)-50:mLocais22(2)+90), 'r','LineWidth',0.2);    
    title('Análise individual das M-waves pré-CVIM, 2ª');
    
    %salva 2 m-waves pós protocolo - segunda
   saveas(fig23, '2 m-waves pós protocolo - segunda.fig')
    
fig24 = figure('units','normalized','outerposition',[0 0 1 1]);
addToolbarExplorationButtons(fig24);
 plot(graficoMwave2(mLocais2(3)-50:mLocais2(3)+90), 'r','LineWidth',0.2); 
    title('Análise individual das M-waves pré-CVIM, 3ª');   
    
    %salva 2 m-waves pós protocolo - terceira
   saveas(fig24, '2 m-waves pós protocolo - terceira.fig')
        
%     uiwait(msgbox({'A análise das m-waves é avaliador-dependente' ''...
%     'Utilize o "Data cursor" para a análise' '' 'Botão com uma caixinha amarela e uma curva passando por ela'}...
%     ,'Análise m-waves','modal'));


% tLocais2135 = [tLocais2(1), tLocais2(3), tLocais2(5)];
% 
% %close all
% %Cálculo e armazenamento das variáveis relevantes
% 
% %   duracaoTq2 = localCrit2' - localOnsetT2(1:3); 
%  duracaoTq2 = localCrit2 - localOnsetT2(1:3); 
%  
%  amplitudeTq2 = tPicos2(1:3) - onsetT2(1:3); 
%  
% %  timetoPeakTq2 = tLocais2(1:3) - localOnsetT2(1:3); 
%  timetoPeakTq2 = tLocais2135(1:3) - localOnsetT2(1:3); 
%  
% %  halfRelaxTimeTq2 = localMeio2 - tLocais2(1:3);  
%  halfRelaxTimeTq2 = localMeio2 - tLocais2135(1:3);  
%  
%  tdf2 = amplitudeTq2 ./ timetoPeakTq2';  
 
 %Caso esteja achando picos a mais, arrumar manualmente os picos
 duracaoTq2 = localCrit2 - localOnsetT2(1:3);
 amplitudeTq2 = tPicos2(1:3) - onsetT2(1:3)';
 timetoPeakTq2 = tLocais2(1:3) - localOnsetT2(1:3);
 halfRelaxTimeTq2 = localMeio2 - tLocais2(1:3);  
 tdf2 = amplitudeTq2 ./ timetoPeakTq2;  

%% Algoritmo que seleciona m-waves pós CVIM e as analiza    
    if resposta == 's' 
    
                                                          
 [mPicos3 mLocais3] = findpeaks(taTara(ponto(4)+15000:end),'MinPeakHeight',4,'MinPeakDistance',5000);   
 [tPicos3 tLocais3] = findpeaks(TqTara(ponto(4)+15000:end),'MinPeakProminence',3,'MinPeakDistance',5000);   


 
 graficoMwave3 = taTara(ponto(4)+15000:end);
 
 %Onsets dos twitchs
[mPicos33 mLocais33] = findpeaks(taTara(ponto(4)+15000:end)*-1,'MinPeakHeight',4,'MinPeakDistance',5000);   
    onsetM3 = mPicos33 *-1;
    
%  Acha o onset plotando o gráfico de torque dentro do intervalo ponto(4)+10000
% ao end adicionando +10 ao torque para deixa-lo positivo e achar um ponto
% crítico na curva e ao final subtrai esses +10 para normalizar

graficoT3 = TqTara(ponto(4)+15000:end)+10;
localOnsetT3 = tLocais3;

localOnsetT3 = localOnsetT3';

for iii = 1:9
    
     while abs(graficoT3(localOnsetT3(iii))) - abs(graficoT3 (localOnsetT3(iii) - 10)) > 0.005
localOnsetT3(iii) = localOnsetT3(iii) - 10;     
     end
end
onsetT3 = graficoT3(localOnsetT3(1:9))-10;

% Achar tempo de meio relaxamento
    graficoMeio3 = graficoT3-10;
        graficoMeio3 = graficoT3;   
amplitudeT3 = tPicos3(1:9) - onsetT3(1:9);
    valorInicio3 = tPicos3(1:9);
valorMeio3 = (tPicos3(1:9)-(amplitudeT3/2))+10;   

    valorFim3 = onsetT3;
localMeio3 = tLocais3(1:9);
    
for i3 = 1:9
    while graficoMeio3(localMeio3(i3)) > valorMeio3(i3)   
        localMeio3(i3) = localMeio3(i3) + 1;
    end
end
    
valorMeio3 = valorMeio3-10;  

% Achar o local de termino da contração *NOVO*
localContinua3 = localMeio3;        

    for i33 = 1:3
%  zero = find(~pontoCritico); 
derivada3 = diff(graficoMeio3(localContinua3(i33):localContinua3(i33)+500));
 
    [pontoCritico3(i33), localCrit3(i33)] = min(abs(derivada3));    
    
localCrit3(i33) = localCrit3(i33) + localContinua3(i33);

    end   
pontoCritico3 = graficoMeio3(localCrit3) - 10;


 % Plot para analisar as m-waves pós CVIM
fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
addToolbarExplorationButtons(fig3);
 plot(taTara(ponto(4)+15000:end), 'r');
    hold on   
 plot(TqTara(ponto(4)+15000:end),'k','LineWidth',1);    
 plot(mLocais3(1:3), mPicos3(1:3),'bo');
 plot(tLocais3(1:3), tPicos3(1:3),'b^');    
 plot(localOnsetT3(1:3), onsetT3(1:3) ,'bv');
 plot(mLocais33(1:3), onsetM3(1:3),'bo') 
 plot(localMeio3(1:3), valorMeio3(1:3),'c*'); 
 plot(localCrit3,pontoCritico3,'bv');
    hold off   
    title('Sinal pós CVIM');

  question = {[' Se sim, digite: 0.' sprintf('\n')...
     ' Caso não estejam, qual deve ser a nova amplitude mínima para que esses picos sejam incluidos na função de busca?'...
     sprintf('\n') 'Amplitude atual: 0.8 ' sprintf('\n') 'Ex: 0.7       Ex2: 0.4']};
    titulo = 'Os picos de torque se encontram no lugar certo?';
    dims = [1 70];
    definput = {'0'};
    answer3 = inputdlg(question,titulo,dims,definput);
    resposta3 = str2num(answer3{:});

    %Função para checar se os picos de torque foram selecionados
    %corretamente

   while resposta3 ~= 0
    
 close;      
    
 % Plot para analisar as m-waves pós CVIM
fig3 = figure('units','normalized','outerposition',[0 0 1 1])
addToolbarExplorationButtons(fig3);
 plot(taTara(ponto(4)+15000:end), 'r');
    hold on   
 plot(TqTara(ponto(4)+15000:end),'k','LineWidth',1);    
 plot(mLocais3(1:3), mPicos3(1:3),'bo');
 plot(tLocais3(1:3), tPicos3(1:3),'b^');    
 plot(localOnsetT3(1:3), onsetT3(1:3) ,'bv');
 plot(mLocais33(1:3), onsetM3(1:3),'bo') 
 plot(localMeio3(1:3), valorMeio3(1:3),'c*'); 
 plot(localCrit3,pontoCritico3,'bv');
    hold off   
    title('Sinal pós CVIM');
    
  question = {[' Se sim, digite: 0.' sprintf('\n')...
     ' Caso não estejam, qual deve ser a nova amplitude mínima para que esses picos sejam incluidos na função de busca?'...
     sprintf('\n') 'Amplitude atual: 0.8 ' sprintf('\n') 'Ex: 0.7       Ex2: 0.4']};
    titulo = ' E agora, os picos de torque se encontram no lugar certo?';
    dims = [1 70];
    definput = {'0'};
%   opts.Resize = 'on';
    answer3 = inputdlg(question,titulo,dims,definput) %opts);
    resposta3 = str2num(answer3{:});    
 
%Repete a função de achar os picos de torque relacionados aos twitches  
    
 [mPicos3 mLocais3] = findpeaks(taTara(ponto(4)+15000:end),'MinPeakHeight',5,'MinPeakDistance',1000);   
 [tPicos3 tLocais3] = findpeaks(TqTara(ponto(4)+15000:end),'MinPeakProminence',resposta3,'MinPeakDistance',1000);   

%  while tLocais3 < 3
%     
%      uiwait(msgbox('A nova amplitude não retornou resultados','modal')); 
%      
%      question = {[' Se sim, digite: 0.' sprintf('\n')...
%      ' Caso não estejam, qual deve ser a nova amplitude mínima para que esses picos sejam incluidos na função de busca?'...
%      sprintf('\n') 'Amplitude atual: 0.8 ' sprintf('\n') 'Ex: 0.7       Ex2: 0.4']};
%     titulo = ' E agora, os picos de torque se encontram no lugar certo?';
%     dims = [1 70];
%     definput = {'0'};
% %   opts.Resize = 'on';
%     answer3 = inputdlg(question,titulo,dims,definput) %opts);
%     resposta3 = str2num(answer3{:});    
%      
%     [mPicos3 mLocais3] = findpeaks(taTara(ponto(4)+15000:end),'MinPeakHeight',5,'MinPeakDistance',1000);   
%  [tPicos3 tLocais3] = findpeaks(TqTara(ponto(4)+15000:end),'MinPeakProminence',resposta3,'MinPeakDistance',1000);   
% 
%  end
 
 graficoMwave3 = taTara(ponto(4)+15000:end);
 
 %Onsets dos twitchs
[mPicos33 mLocais33] = findpeaks(taTara(ponto(4)+15000:end)*-1,'MinPeakHeight',3,'MinPeakDistance',1000);   
    onsetM3 = mPicos33 *-1;
    
%  Acha o onset plotando o gráfico de torque dentro do intervalo ponto(4)+10000
% ao end adicionando +10 ao torque para deixa-lo positivo e achar um ponto
% crítico na curva e ao final subtrai esses +10 para normalizar

graficoT3 = TqTara(ponto(4)+15000:end)+10;
localOnsetT3 = tLocais3;
     
for iii = 1:9
     while abs(graficoT3(localOnsetT3(iii))) - abs(graficoT3 (localOnsetT3(iii) - 10)) > 0.0005
localOnsetT3(iii) = localOnsetT3(iii) - 10;     
     end
end
onsetT3 = graficoT3(localOnsetT3(1:9))-10;

% Achar tempo de meio relaxamento
    graficoMeio3 = graficoT3-10;
        graficoMeio3 = graficoT3;   
amplitudeT3 = tPicos3(1:9) - onsetT3(1:9);
    valorInicio3 = tPicos3(1:9);
valorMeio3 = (tPicos3(1:9)-(amplitudeT3/2))+10;   

    valorFim3 = onsetT3;
localMeio3 = tLocais3(1:9);
    
for i3 = 1:9
    while graficoMeio3(localMeio3(i3)) > valorMeio3(i3)   
        localMeio3(i3) = localMeio3(i3) + 1;
    end
end
    
valorMeio3 = valorMeio3-10;  

% Achar o local de termino da contração *NOVO*
localContinua3 = localMeio3;        

    for i33 = 1:3
%  zero = find(~pontoCritico); 
derivada3 = diff(graficoMeio3(localContinua3(i33):localContinua3(i33)+500));
 
    [pontoCritico3(i33), localCrit3(i33)] = min(abs(derivada3));    
    
localCrit3(i33) = localCrit3(i33) + localContinua3(i33);

    end   
pontoCritico3 = graficoMeio3(localCrit3) - 10;  
    
   end   
   
    %salva 3 twitches pós CVIM
   saveas(fig3, '3 twitches pós CVIM.fig')
   
%     uiwait(msgbox('Será feita a análise dos 3 primeiros Twitchs pós CVIM','Twitchs','modal'));
 
%plot das 3 m-waves pós CVIM
     fig31= figure('units','normalized','outerposition',[0 0 1 1]);
     addToolbarExplorationButtons(fig31);
subplot(3,1,1)
    plot(graficoMwave3(mLocais3(1)-50:mLocais3(1)+90), 'r','LineWidth',0.2);
    title('Análise individual das M-waves pós-CVIM, 1ª');  
subplot(3,1,2)        
    plot(graficoMwave3(mLocais33(2)-50:mLocais33(2)+90), 'r','LineWidth',0.2);    
    title('Análise individual das M-waves pós-CVIM, 2ª');
subplot(3,1,3) 
    plot(graficoMwave3(mLocais3(3)-50:mLocais3(3)+90), 'r','LineWidth',0.2); 
    title('Análise individual das M-waves pós-CVIM, 3ª');     
 
    %salva 3 twitches pós CVIM
   saveas(fig31, '3 twitches pós CVIM, os 3.fig')
    
%     
%  uiwait(msgbox({'A análise das m-waves é avaliador-dependente' ''...
%     'Utilize o "Data cursor" para a análise' '' 'Botão com uma caixinha amarela e uma curva passando por ela'}...
%     ,'Análise m-waves','modal')); 


%plot da primeira m-wave pós CVIM
fig32 = figure('units','normalized','outerposition',[0 0 1 1]);
addToolbarExplorationButtons(fig32);
 plot(graficoMwave3(mLocais3(1)-50:mLocais3(1)+90), 'r','LineWidth',0.2);
    title('Análise individual das M-waves pós-CVIM, 1ª');  
    
    %salva 3 m-waves pós CVIM - primeira
   saveas(fig32, '3 m-waves pós CVIM - primeira.fig')
   
%plot da segunda m-wave pós CVIM   
fig33 = figure('units','normalized','outerposition',[0 0 1 1]);
addToolbarExplorationButtons(fig33);
 plot(graficoMwave3(mLocais33(2)-50:mLocais33(2)+90), 'r','LineWidth',0.2);    
    title('Análise individual das M-waves pós-CVIM, 2ª');
    
     %salva 3 m-waves pós CVIM - segunda
   saveas(fig33, '3 m-waves pós CVIM - segunda.fig')

%plot da terceira m-wave pós CVIM   
fig34 = figure('units','normalized','outerposition',[0 0 1 1]);
addToolbarExplorationButtons(fig34);
 plot(graficoMwave3(mLocais3(3)-50:mLocais3(3)+90), 'r','LineWidth',0.2); 
    title('Análise individual das M-waves pós-CVIM, 3ª');  
    
    %salva 3 m-waves pós CVIM - terceira
   saveas(fig34, '3 m-waves pós CVIM - terceira.fig')
    
%     uiwait(msgbox({'A análise das m-waves é avaliador-dependente' ''...
%     'Utilize o "Data cursor" para a análise' '' 'Botão com uma caixinha amarela e uma curva passando por ela'}...
%     ,'Análise m-waves','modal'));,

    
    %close all  
%Cálculo e armazenamento das variáveis relevantes
 duracaoTq3 = localCrit3' - localOnsetT3(1:3);
 amplitudeTq3 = tPicos3(1:3) - onsetT3(1:3);
 timetoPeakTq3 = tLocais3(1:3) - localOnsetT3(1:3);
 halfRelaxTimeTq3 = localMeio3(1:3) - tLocais3(1:3);  
 tdf3 = amplitudeTq3 ./ timetoPeakTq3;   
 
    end  
    
    %Achar o local de termino da contração *JEITO DOIDO DO JAPA*

%     localContinua3 = localMeio3;        
% for i33 = 1:9
%  
%     media3(i33) = mean(graficoMeio3(localContinua3(i33)+400:localContinua3(i33)+700));    
%     Desvio3(i33) = 3*std(graficoMeio3(localContinua3(i33)+400:localContinua3(i33)+700));        
%     mediaDesvio3(i33) = media3(i33)+Desvio3(i33);
% 
%     if graficoMeio3(localContinua3(i33)+500) >= valorFim3(i33)         
%                localContinua3(i33) = localContinua3(i33)+500;
%               while graficoMeio3(localContinua3(i33)) <= mediaDesvio3(i33) 
%               localContinua3(i33) = localContinua3(i33) -1;
%               end
%     else 
%               while graficoMeio3(localContinua3(i33)) > valorFim3(i33)  
%     localContinua3(i33) = localContinua3(i33) + 1;
%               end      
%     end
% end
%% Divide sinal em duas funções  

%acha os picos máximos da função
[picos,localmax] = findpeaks(torqueCortadoPos(:),'MinPeakDistance',1200);

%acha os picos mínimos da função
[picosmin,localmin] = findpeaks(torqueCortadoPos(:)*-1,'MinPeakDistance',1200);

%o picos minimos voltam ao original
picosminOrig = picosmin*-1;

% Deixa valores positivos
% picos = picos + menorVal;
% picosminOrig = picosmin*-1 + menorVal;
 

%% Root mean square EMG

%Eletro-estimulações
% EMGsquare1 = taTara(ponto(2):ponto(3)).^2;
% mediaEMG1 = mean(EMGsquare1);
%     rmsEMG1 = sqrt(mediaEMG1); 
%  
% fprintf('A rms da eletro-miografia das 250 contrações é: %.2f microV \n', rmsEMG1)    
%% Root mean square CVM platô

    if resposta == 's'

EMGsquare2 = taTara(platox1-1001:platox1+1000).^2; 
mediaEMG2 = mean(EMGsquare2); 
    rmsEMG2 = sqrt(mediaEMG2); 

fprintf('A rms da eletro-miografia do platô de torque na CVIM é: %f mV \n', rmsEMG2)    
    end

%% Calcula área funções

areasup = trapz(picos);
% fprintf('A área superior é: %.2f \n',areasup)

areainf = trapz(picosminOrig);
% fprintf('A área inferior é: %.2f \n', areainf)

difarea = areasup - areainf;
fprintf('A área contida no vetor de torque (Impulso) é: %.2f N.s \n', difarea)

    if resposta == 's'
fprintf('A média de torque do platô é: %.2f N \n', mediaPlato)
    end

%% Pico a pico
    %Achar o local dos onsets
% localOnset1 = max(1:localmax(1));
% localOnset = localmax(2:end) - 1000; 
% 
% %Achar o valor desses onsets (todos devem ser positivos)
% onset1 = torqueCortado(localOnset1);
% onset = torqueCortadoPos(localOnset);
 
if length(picos) <= length(picosmin)
    compri = length(picos);

elseif length(picosmin) <= length(picos)
    compri = length(picosmin);
end 

%calcula a amplitude de cada contração (250 contrações)
amplitude = picos(1:compri) - picosminOrig(1:compri);

%% Pico a pico em setores com incrementos

  [picossetorial1] = calculaPicos(picos,picosminOrig,1);
  [picossetorial5] = calculaPicos(picos,picosminOrig,5);
  [picossetorial10] = calculaPicos(picos,picosminOrig,10);
  [picossetorial25] = calculaPicos(picos,picosminOrig,25);
  [picossetorial50] = calculaPicos(picos,picosminOrig,50);
  
%% Plot dos gráficos
 fig4 = figure('units','normalized','outerposition',[0 0 1 1]);
 addToolbarExplorationButtons(fig4);
%Plota o sinal original ja descontando o peso da perna (começa em zero)
subplot(3,1,1)
 
    plot(taTara, 'r','LineWidth',0.2)
    hold on
    plot(TqTara,'k','LineWidth',1)
    title('Sinal de torque e EMG (com baseline)')

%Plota duas linhas vermelhas tracejadas no local de corte da
% função ginput
y1=get(gca,'ylim');
    plot([ponto(2) ponto(2)],y1,'r-.')
    plot([ponto(3) ponto(3)],y1,'r-.')
  
subplot(3,1,2)
%Plota valores máximos e mínimos já positivados
    
     plot(picossetorial1,'k','LineWidth',1)
hold on
    plot(picos)
    plot(picosminOrig,'g');  
    title('Valores máx e min positivos de torque');
    
%Plota o trabalho (W) em janelas de 25, convertido para porcentagem
subplot(3,1,3)

picosMax = max(picossetorial25);
picosplot = (picossetorial25 / picosMax)*100;

    plot(picosplot)
% plot(areasetorial25)
    title('Valores de torque N (Janelas de 25)');
    xlabel('Tempo')
    ylabel('Porcentagem %')
    grid on
hold off

%salva Gráficos
   saveas(fig4, 'Gráficos.fig')

%% Plot dos valores de pico a pico N.m
fig5 = figure('units','normalized','outerposition',[0 0 1 1]);
addToolbarExplorationButtons(fig5);

subplot(5,1,1)
    plot(picossetorial1,'k','LineWidth',1)
hold on
    plot(picos,'--')
    plot(picosminOrig,'g--');  
    title('Os 250 valores de torque N ')
subplot(5,1,2)
    plot(picossetorial5,'b')
    title('50 valores')
subplot(5,1,3)
    plot(picossetorial10,'m')
    title('25 valores')
subplot(5,1,4)
    plot(picossetorial25,'r')
    title('10 valores')
subplot(5,1,5)
    plot(picossetorial50,'c')
    title('5 valores')
hold off  

%salva Gráficos em janelas
   saveas(fig5, 'Gráficos em janelas.fig')

%% Trabalho
%Porque dá um valor tão alto?
%Dividir por 10000 pra ser mais coerente
trabalho = trapz(torqueCortadoPos);
tbb = round(trabalho)/10000;

fprintf('O impulso realizado dentro das 250 contrações é: %.2f N.s \n', tbb)


%% Calcula média da m-wave durante platô da CVIM

%plot(taTara(ponto(4):ponto(4)+22500),'k','LineWidth',1);

periodoMwave = (platox1-500:platox1+500);
valoresMwave = taTara(platox1-500:platox1+500);
mediaMwaveCVIM = rms(valoresMwave);

fprintf('A média da eletromiografia durante o platô da CVIM (utilizando função rms) é: %.2f mV \n', mediaMwaveCVIM);

% 
% periodoMwave = (platox-500:platox+500);
% valoresMwave = ta(platox-500:platox+500);
% mediaMwaveCVIM = rms(valoresMwave);
% 
% fprintf('A média da eletromiografia durante o platô da CVIM (utilizando função rms) é: %.2f mV \n', mediaMwaveCVIM);


%% Plot das 9 m-waves entrepostas

fig6 = figure('units','normalized','outerposition',[0 0 1 1]);
plot(taTara(mLocais(3)-86:mLocais(3)+90), 'b','LineWidth',1.5);
	hold on 
plot(taTara(mLocais(2)-50:mLocais(2)+90), 'b','LineWidth',1.5,'HandleVisibility','off'); 
plot(taTara(mLocais(1)-50:mLocais(1)+90), 'b','LineWidth',1.5,'HandleVisibility','off');
	
	plot(graficoMwave2(mLocais2(1)-84:mLocais2(1)+90), 'g','LineWidth',1.5);
	plot(graficoMwave2(mLocais22(2)-52:mLocais22(2)+90), 'g','LineWidth',1.5,'HandleVisibility','off');
	plot(graficoMwave2(mLocais2(3)-50:mLocais2(3)+90), 'g','LineWidth',1.5,'HandleVisibility','off');
	
		plot(graficoMwave3(mLocais3(1)-80:mLocais3(1)+90), 'r','LineWidth',1.5)
		plot(graficoMwave3(mLocais3(2)-82:mLocais3(2)+90), 'r','LineWidth',1.5,'HandleVisibility','off')
		plot(graficoMwave3(mLocais3(3)-50:mLocais3(3)+90), 'r','LineWidth',1.5,'HandleVisibility','off')

legend('Pré-protocolo','Pós-protocolo','Pós-CVIM');

%% 
1+1 == 2;

%% Salva arquivos
if resposta == 's'
save([filename(1:end-4),'.mat'],'Torquefilt','ponto','torqueCortado',...
    'picos','picosmin','tbb','amplitude','picossetorial5',...
    'picossetorial10','picossetorial25','picossetorial50','duracaoTq1',...
    'amplitudeTq1','timetoPeakTq1' , 'halfRelaxTimeTq1', 'tdf1', 'mediaPlato','rmsEMG2',...
    'duracaoTq2', 'amplitudeTq2', 'timetoPeakTq2', 'halfRelaxTimeTq2', 'tdf2',...
    'duracaoTq3', 'amplitudeTq3', 'timetoPeakTq3', 'halfRelaxTimeTq3', 'tdf3');

% saveas(fig1,[filename(1:end-4),'_1.png']);
% 
% saveas(fig2,[filename(1:end-4),'_2.png']);
% 
% saveas(fig3,[filename(1:end-4),'_3.png']);
% 
% saveas(fig4,[filename(1:end-4),'_4.png']);
% 
% saveas(fig5,[filename(1:end-4),'_5.png']);

gg = load([filename(1:end-4),'.mat']);

nomes = {'ponto','picos','picosmin','Trabalho','pico a pico (torque)', 'Torque setores de 5',...
    'Torque setores de 10','Torque setores de 25','Torque setores de 50'...
    'Duracao torque pré','Amplitude torque pré','Tempo pro pico pré','Tempo de meio relax pré',...
    'TDF pros picos pré','Média platô CVIM','RMS da EMG no platô CVIM', ...
    'Duracao torque pós','Amplitude torque pós','Tempo pro pico pós','Tempo de meio relax pós',...
    'TDF pros picos pós', 'Duracao torque pós2','Amplitude torque pós2','Tempo pro pico pós2',...
    'Tempo de meio relax pós2','TDF pros picos pós2' }; 

nomes{2,1} = ponto;
nomes{2,2} = picos;
nomes{2,3} = picosmin;
nomes{2,4} = tbb;
nomes{2,5} = amplitude;
nomes{2,6} = picossetorial5';
nomes{2,7} = picossetorial10';
nomes{2,8} = picossetorial25';
nomes{2,9} = picossetorial50';
nomes{2,10} = duracaoTq1;
nomes{2,11} = amplitudeTq1;
nomes{2,12} = timetoPeakTq1;
nomes{2,13} = halfRelaxTimeTq1;
nomes{2,14} = tdf1;
nomes{2,15} = mediaPlato;
nomes{2,16} = rmsEMG2;
nomes{2,17} = duracaoTq2;
nomes{2,18} = amplitudeTq2;
nomes{2,19} = timetoPeakTq2;
nomes{2,20} = halfRelaxTimeTq2;
nomes{2,21} = tdf2;
nomes{2,22} = duracaoTq3;
nomes{2,23} = amplitudeTq3;
nomes{2,24} = timetoPeakTq3;
nomes{2,25} = halfRelaxTimeTq3;
nomes{2,26} = tdf3;

xlswrite([filename(1:end-4),'.xls'],nomes,1,'A1')
xlswrite([filename(1:end-4),'.xls'],nomes{2,1},1,'A2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,2},1,'B2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,3},1,'C2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,4},1,'D2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,5},1,'E2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,6},1,'F2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,7},1,'G2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,8},1,'H2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,9},1,'I2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,10},1,'J2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,11},1,'K2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,12},1,'L2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,13},1,'M2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,14},1,'N2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,15},1,'O2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,16},1,'P2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,17},1,'Q2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,18},1,'R2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,19},1,'S2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,20},1,'T2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,21},1,'U2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,22},1,'V2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,23},1,'W2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,24},1,'X2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,25},1,'Y2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,26},1,'Z2')

%Tentativa de largura da coluna
% hExcel = actxserver('Excel.Application')
% hWorkbook = hExcel.Workbooks.Open(filenamee)
% hWorksheet = hWorkbook.Sheets.Item(1)
% hWorksheet.Columns.Item(1).columnWidth = 100; %first column
% hWorksheet.Columns.Item(2).columnWidth = 100; %second column

%%
elseif resposta == 'n'
    
save([filename(1:end-4),'.mat'],'Torquefilt','ponto','torqueCortado',...
    'picos','picosmin','tbb','amplitude','picossetorial5',...
    'picossetorial10','picossetorial25','picossetorial50','duracaoTq1',...
    'amplitudeTq1','timetoPeakTq1' , 'halfRelaxTimeTq1', 'tdf1',...
    'duracaoTq2', 'amplitudeTq2','timetoPeakTq2' , 'halfRelaxTimeTq2', 'tdf2');

% saveas(fig1,[filename(1:end-4),'_1.png']);
% 
% saveas(fig2,[filename(1:end-4),'_2.png']);
% 
% saveas(fig4,[filename(1:end-4),'_4.png']);
% 
% saveas(fig5,[filename(1:end-4),'_5.png']);

gg = load([filename(1:end-4),'.mat']);

nomes = {'ponto','picos','picosmin','Trabalho','pico a pico (torque)', 'Torque setores de 5',...
    'Torque setores de 10','Torque setores de 25','Torque setores de 50'...
    'Duracao torque pré','Amplitude torque pré','Tempo pro pico pré','Tempo de meio relax pré',...
    'TDF pros picos pré', 'Duracao torque pós','Amplitude torque pós','Tempo pro pico pós','Tempo de meio relax pós',...
    'TDF pros picos pós'}; 

nomes{2,1} = ponto;
nomes{2,2} = picos;
nomes{2,3} = picosmin;
nomes{2,4} = tbb;
nomes{2,5} = amplitude;
nomes{2,6} = picossetorial5';
nomes{2,7} = picossetorial10';
nomes{2,8} = picossetorial25';
nomes{2,9} = picossetorial50';
nomes{2,10} = duracaoTq1;
nomes{2,11} = amplitudeTq1;
nomes{2,12} = timetoPeakTq1;
nomes{2,13} = halfRelaxTimeTq1;
nomes{2,14} = tdf1;
nomes{2,15} = duracaoTq2;
nomes{2,16} = amplitudeTq2;
nomes{2,17} = timetoPeakTq2;
nomes{2,18} = halfRelaxTimeTq2;
nomes{2,19} = tdf2;


xlswrite([filename(1:end-4),'.xls'],nomes,1,'A1')
xlswrite([filename(1:end-4),'.xls'],nomes{2,1},1,'A2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,2},1,'B2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,3},1,'C2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,4},1,'D2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,5},1,'E2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,6},1,'F2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,7},1,'G2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,8},1,'H2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,9},1,'I2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,10},1,'J2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,11},1,'K2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,12},1,'L2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,13},1,'M2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,14},1,'N2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,15},1,'O2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,16},1,'P2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,17},1,'Q2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,18},1,'R2')
xlswrite([filename(1:end-4),'.xls'],nomes{2,19},1,'S2')

end
end

%%Salva arquivos no caso de CVIM existente e falta de twitches pós CVIM

% save([filename(1:end-4),'.mat'],'Torquefilt','ponto','torqueCortado',...
%     'picos','picosmin','tbb','amplitude','picossetorial5',...
%     'picossetorial10','picossetorial25','picossetorial50','duracaoTq1',...
%     'amplitudeTq1','timetoPeakTq1' , 'halfRelaxTimeTq1', 'tdf1', 'mediaPlato','rmsEMG2',...
%     'duracaoTq2', 'amplitudeTq2', 'timetoPeakTq2', 'halfRelaxTimeTq2', 'tdf2');
% 
% gg = load([filename(1:end-4),'.mat']);
% 
% nomes = {'ponto','picos','picosmin','Trabalho','pico a pico (torque)', 'Torque setores de 5',...
%     'Torque setores de 10','Torque setores de 25','Torque setores de 50'...
%     'Duracao torque pré','Amplitude torque pré','Tempo pro pico pré','Tempo de meio relax pré',...
%     'TDF pros picos pré','Média platô CVIM','RMS da EMG no platô CVIM', ...
%     'Duracao torque pós','Amplitude torque pós','Tempo pro pico pós','Tempo de meio relax pós',...
%     'TDF pros picos pós'}; 
% 
% nomes{2,1} = ponto;
% nomes{2,2} = picos;
% nomes{2,3} = picosmin;
% nomes{2,4} = tbb;
% nomes{2,5} = amplitude;
% nomes{2,6} = picossetorial5';
% nomes{2,7} = picossetorial10';
% nomes{2,8} = picossetorial25';
% nomes{2,9} = picossetorial50';
% nomes{2,10} = duracaoTq1;
% nomes{2,11} = amplitudeTq1;
% nomes{2,12} = timetoPeakTq1;
% nomes{2,13} = halfRelaxTimeTq1;
% nomes{2,14} = tdf1;
% nomes{2,15} = mediaPlato;
% nomes{2,16} = rmsEMG2;
% nomes{2,17} = duracaoTq2;
% nomes{2,18} = amplitudeTq2;
% nomes{2,19} = timetoPeakTq2;
% nomes{2,20} = halfRelaxTimeTq2;
% nomes{2,21} = tdf2;
% 
% xlswrite([filename(1:end-4),'.xls'],nomes,1,'A1')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,1},1,'A2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,2},1,'B2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,3},1,'C2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,4},1,'D2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,5},1,'E2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,6},1,'F2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,7},1,'G2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,8},1,'H2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,9},1,'I2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,10},1,'J2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,11},1,'K2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,12},1,'L2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,13},1,'M2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,14},1,'N2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,15},1,'O2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,16},1,'P2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,17},1,'Q2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,18},1,'R2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,19},1,'S2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,20},1,'T2')
% xlswrite([filename(1:end-4),'.xls'],nomes{2,21},1,'U2')


%% Função que calcula amplitudes em diferentes incrementos
function [picossetorial] = calculaPicos(picos,picosminOrig,inc)

%incrementos 'inc' definidos la em cima
inicioo = 1;
picossetorial1=[];

while inicioo <= length(picos)-(inc)
    
    picossetorial1 = [picossetorial1, mean(picos(inicioo:inicioo+inc))];
    
    inicioo = inicioo+inc;
    
end

inicioo2 = 1;
picossetorial2 = [];

while inicioo2 <= length(picosminOrig)-(inc)
    
    picossetorial2 = [picossetorial2, mean(picosminOrig(inicioo2:inicioo2+inc))];
    
    inicioo2 = inicioo2+inc;
    
end
    
if length(picossetorial1) <= length(picossetorial2)
    limitee = length(picossetorial1);

elseif length(picossetorial2) <= length(picossetorial1)
    limitee = length(picossetorial2);
end

picossetorial = picossetorial1(1:limitee) - picossetorial2(1:limitee);

end
%% Função que calcula a area(Trabalho W) dividinto o torque em setores(inc)
function [areasetorial] = calculaArea(torqueCortadoPos, inc)

inicio1 = 1;
areasetorial1=[];
%inc = 25;

while inicio1 < length(torqueCortadoPos)-(inc)
   
  areasetorial1 = [areasetorial1 , trapz(torqueCortadoPos(inicio1:(inicio1+inc)))];

  inicio1 = inicio1+inc;

end

areasetorial = areasetorial1;
end



function [x,y] = ginput_zoom(w)
%
% [x,y] = ginput_zoom()
%
%  You can zoom and pan before selecting each point,
%  press 'space' when you are done to start the selection
%
%  If you made a mistake, you can delete latest points pressing 'backspace'
% 
%  Once you are finished, press 'enter' to exit
%
%  Date: 27/01/2017
%  Authors: Martin Sanz Sabater (University of Valencia)
%           Javier Garcia Monreal (University of Valencia)
%

x=[];y=[];
h=[];
w=1;
while w~=0; 
  w = waitforbuttonpress;
  while w==0
    w = waitforbuttonpress;
  end
  cfg = gcf();
  ch = double(get(cfg, 'CurrentCharacter'));
  if ch == 13 % ENTER button
    break;
  end
  if ch == 8 % Backspace button
    if isempty(x) == 0
        x = x(1:end-1);
        y = y(1:end-1);
        delete(h(end));
        h = h(1:end-1);
        continue;
    end
  end
  
  [a,b]=ginput(1);
  x = [x;a];
  y = [y;b];
  hold on; h = [h;plot(x,y,'g+')]; hold off;
end
%%
 function  [X,Y,BUTTON,SCALEMAT] = ginput2(varargin)
% %GINPUT2   Graphical input from mouse with zoom, pan, plot and scaling.
% %   
% %   SYNTAX:
% %                        XY = ginput2;
% %                        XY = ginput2(DoScale);          true or false
% %                        XY = ginput2(...,PlotOpt);      '.r' for example
% %                        XY = ginput2(...,'KeepZoom');   vs. 'UnZoom'
% %                        XY = ginput2(N,...);
% %                        XY = ginput2(...);
% %                     [X,Y] = ginput2(...);
% %              [X,Y,BUTTON] = ginput2(...);
% %     [X,Y,BUTTON,SCALEMAT] = ginput2(...);
% %
% %   INPUT:
% %     DoScale    - Single logical specifying whether the IMAGE should be
% %                  interactively scaled (georeferenced), or it can be the
% %                  2x4 SCALEMAT matrix for automatically scaling.
% %                  DEFAULT: false (do not scales/georeferences)
% %     PlotOpt    - String and/or parameter/value pairs specifying the drawn
% %                  points optional inputs (see PLOT for details). 
% %                  DEFAULT: 'none' (do not plot any point)
% %     'KeepZoom' - When finishing selection by default the zoom is
% %                  restored. By using this option this feature is ignored.
% %                  DEFAULT: 'UnZoom' (restores original axis limits)
% %     N          - Number of points to be selected. One of 0,1,2,...,Inf
% %                  DEFAULT: Inf (selects until ENTER or ESCAPE is pressed)
% %
% %   OUTPUT:
% %     XY        - [X(:) Y(:)] axis coordinate(s).
% %     X         - X-coordinate(s).
% %     Y         - Y-coordinate(s).
% %     BUTTON    - Last pressed button.
% %     SCALEMAT  - 2x4 matrix specifying the coordinates of two different
% %                 points (1) and (2) in the Image coordinates (pixels) and
% %                 the User coordinates (data):
% %                                                Point 1     Point 2
% %                   Image coord (pixels):     [ (I1x,I1y)   (I2x,I2y) 
% %                   User  coord (data)  :       (U1x,U1y)   (U2x,U2y) ]
% %                 to be use for scaling/georeferencing.
% %
% %   DESCRIPTION:
% %     This program uses MATLAB's GINPUT function to get the coordinates
% %     of a mouse-selected point in the current figure (see GINPUT for
% %     details), but with five major improvements:
% %                  1. ZOOMING  (left click)
% %                  2. PANNING  (dragging mouse)
% %                  3. DELETING (last selected point)
% %                  4. PLOTING  (temporarily the selected points)
% %                  5. SCALING or GEOREFERENCE IMAGES.
% %     The differences are:
% %      a) Obviously, the SCALEOPT and PlotOpt optional arguments.
% %      b) When click is made outside the axes, it is ignored.
% %      c) When LEFT-click, ZOOM-IN is performed right into the selected
% %         point (PANNING).
% %      d) When RIGHT-click, the point is selected (normal).
% %      e) When DOUBLE-click, ZOOM-OUT is done.
% %      f) When MIDDLE-click, ZOOM-RESET is done (see ZOOM for details).
% %      g) When dragging while pressed left-click PAN is done (until the
% %         button is released).
% %      h) When pressed any KEY follows the next rules: 
% %          A) If ENTER is pressed, the selection is terminated. If no point
% %             was already selected, the outputs are empty's.
% %          B) If BACKSPACE key is pressed, the last selected point is
% %             deleted and the selection continues.
% %          C) If SPACEBAR the mouse current position or NANs coordinates
% %             are saved, depending whether the mouse was inside or outside
% %             any of the current figure axes, respectively. In this latter
% %             case, the selection is NOT counted as one of the N points.
% %             Besides, when drawing the color is changed. Then, the outputs
% %             may not be of length N.
% %
% %   NOTE:
% %     * Optional inputs use its DEFAULT value when not given or [].
% %     * Optional outputs may or not be called.
% %     * String inputs may be shortened, as long as they are unambiguous.
% %       Case is ignored.
% %     * The function can be used for interactively digitalize/vectorize
% %       RASTER images with:
% %       >> ginput(true)
% %     * The function can be used only as a georeference function with 
% %       >> ginput2(0,true)
% %     * The scale/georeference only works when the current axes has an
% %       IMAGE type children (see Image for details). 
% %     * The x and y data from axes and image are changed when scale/
% %       georeference is used.
% %     * The drawn points are deleted from the graphics once the selection
% %       is finished. 
% %     * The priority of the inputs are: N, then SCALEOPT and finally
% %       PlotOpt. If the first (integer) is missing, the next is taken into
% %       account (logical or 2x4 matrix) and so on.
% %
% %   EXAMPLE:
% %     % Selects until ENTER is pressed:
% %         xy = ginput2;
% %     % Selects 5 points:
% %         [x,y] = ginput2(5);
% %     % Gets pressed button:
% %         [x,y,button] = ginput2(1);
% %     % Scales image and select 4 points temporarily coloring them in
% %       black. Besides to not ZOOM OUT at the end:
% %         imagesc(peaks(40))
% %         [x,y,button,scalemat] = ginput2(4,true,'k*','KeepZoom');
% %         hold on, plot(x,y,'or'), hold off
% %
% %   SEE ALSO:
% %     GINPUT, PLOT.
% %
% %
% %   ---
% %   MFILE:   ginput2.m
% %   VERSION: 3.1 (Nov 12, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
% %   MATLAB:  7.7.0.471 (R2008b)
% %   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
% %   CONTACT: nubeobscura@hotmail.com
% 
% %   REVISIONS:
% %   1.0      Released. (Jul 09, 2008)
% %   2.0      Changed default YESERROR value and fixed a bug there. Changed
% %            behavior when N==1. Fixed bug with zoom out. Changed default
% %            selection keys. Changed default selection click mouse: from
% %            left one to the right one. (Jun 08, 2009)
% %   2.1      Fixed bugs related with points deletion. Added new 'KeepZoom'
% %            feature. (Aug 20, 2009)
% %   3.0      Now it PANs when dragging. Updated help. (Nov 05, 2009)
% %   3.1      Now returns when N==1 and pressed not predefined KEYS or one
% %            of DELECTION or RETURN buttons. (Nov 12, 2009)
% 
% %   DISCLAIMER:
% %   ginput2.m is provided "as is" without warranty of any kind, under the
% %   revised BSD license.
% 
% %   Copyright (c) 2008-2009 Carlos Adrian Vargas Aguilera
% 
% 
% % INPUTS CHECK-IN
% % -------------------------------------------------------------------------
% 
% % PARAMETERS
% % Defaults:
% X        = [];
% Y        = [];
% BUTTON   = [];
% SCALEMAT = [];
% N        = Inf;
% DoScale  = false;
% PlotOpt  = {'none'};
% UnZoom   = 'UnZoom';
% % Constants KEYs (on my personal keyboard):
% DOUBLECLICK    =   0;
% LEFTCLICK      =   1;
% MIDDLECLICK    =   2;
% RIGHTCLICK     =   3;
% BACKSPACE      =   8;
% ESCAPE         =  27;
% LEFTARROW      =  28;
% RIGHTARROW     =  29;
% UPARROW        =  30;
% DOWNARROW      =  31;
% SPACEBAR       =  32;
% DELETE         = 127;
% ASCII          = [ ...
%                     33:64  ...  UP-KEYS
%                     65:90  ...  UP-LETTERS
%                     91:96  ... LOW-KEYS
%                     97:122 ... LOW-LETTERS
%                    123:126 ... LOW-KEY
%                    161:255 ...     FOREING
%                    ];
% % Functionality:
% % NOTE: I left all this KEYs because the user may use this special case for
% % other special purposes outside this function.
% % % First version default:
% % % SELECTS   = [LEFTCLICK ASCII ESCAPE LEFTARROW RIGHTARROW ...
% % %              UPARROW DOWNARROW SPACEBAR DELETE];
% % % ZOOMIN    = RIGHTCLICK; 
% SELECTS   = [RIGHTCLICK SPACEBAR]; % Selection  buttons
% DELETES   = BACKSPACE;             % Deletion   buttons
% FINISHES  = [];                    % Finishes   buttons
% ZOOMIN    = LEFTCLICK;             % ZOOM(2)    buttons
% ZOOMRESET = MIDDLECLICK;           % ZOOM RESET buttons
% ZOOMOUT   = DOUBLECLICK;           % ZOOM OUT   buttons
% % Other parameters
% secpause  = 0.3;    % Seconds to wait for double-click response.
% YESERROR  = false;  % If there is an error with GINPUT, it tells to display 
%                     % an ERROR or a WARNING message.
%                     
% % Checks number of inputs:
% if nargout>4
%  error('CVARGAS:ginput2:tooManyOutputs',...
%   'At most 4 outputs are allowed.')
% end
% 
% % Checks N:
% if ~isempty(varargin) && ~isempty(varargin{1}) && ...
%   isfloat(varargin{1}) 
%  N           = round(abs(varargin{1}(1))); % Forced unique, positive 
%  varargin(1) = [];                         % integer.
% end
% 
% % Checks DoScale:
% if ~isempty(varargin) && ~isempty(varargin{1}) && ...
%    ((islogical(varargin{1})) || (ndims(varargin{1})==2 && ...
%     all(size(varargin{1})==[2 4])))
%  DoScale     = varargin{1};
%  varargin(1) = [];
% end
% 
% % Checks UnZoom:
% if ~isempty(varargin) 
%  if ~isempty(varargin{1}) && ischar(varargin{1})
%   if     strncmpi(varargin(1),'UnZoom'  ,max(length(varargin{1}),2))
%    UnZoom = 'UnZoom';
%    varargin(1) = [];
%   elseif strncmpi(varargin(1),'KeepZoom',max(length(varargin{1}),2))
%    UnZoom = 'KeepZoom';
%    varargin(1) = [];
%   end
%  elseif (length(varargin)>1) && ~isempty(varargin{end}) && ...
%      ischar(varargin{end})
%   if     strncmpi(varargin(end),'UnZoom'  ,max(length(varargin{1}),2))
%    UnZoom = 'UnZoom';
%    varargin(end) = [];
%   elseif strncmpi(varargin(end),'KeepZoom',max(length(varargin{1}),2))
%    UnZoom = 'KeepZoom';
%    varargin(end) = [];
%   end
%  end
% end
% 
% % Checks PlotOpt:
% if ~isempty(varargin) && ~isempty(varargin{1})
%  PlotOpt = varargin;
% end
% clear varargin
% 
% % Checks DoScale:
% if ~islogical(DoScale)
%  SCALEMAT = DoScale;
%  DoScale = true;
% end
% 
% % SCALES/GEOREFERENCE?:
% if DoScale
%  method = 'linear';
%  extrap = 'extrap';
%  ha     = gca;
%  hi     = findobj(get(ha,'Children'),'Type','image');
%  axes(ha)
%  if ~isempty(hi)
%   hi    = hi(1);
%   xlim  = get(ha,'XLim');
%   ylim  = get(ha,'YLim');
%   zlim  = get(ha,'ZLim');
%   z     = repmat(max(zlim),1,5);
%   xdata = get(hi,'XData');
%   ydata = get(hi,'YData');
%   if isempty(SCALEMAT) % interactively
%    I1x = round(min(xdata)); I2x = round(max(xdata));
%    I1y = round(min(ydata)); I2y = round(max(ydata));
%    % Default (equal):
%    U1x = I1x; U2x = I2x;
%    U1y = I1y; U2y = I2y;
%    hgeo     = [];
%    dlgTitle = 'Georeference image';
%    lineNo   = 1;
%   
%    while true
%     % Selects first corner:
%     theans = ...
%           questdlg('Select the first corner (1 of 2):',dlgTitle,'OK','OK');
%     if ~strcmp(theans,'OK'), return, end
%     pause(secpause)
%     
%     [I1x,I1y] = ginput2(1,false,'none','UnZoom');
%     I1x       = round(I1x);
%     I1y       = round(I1y);
%     if ~ishandle(ha), return, end
%     if (ha==gca) && ~isempty(I1x) && ~isnan(I1x)
%      axis(ha,[xlim ylim])
%      hgeo(1) = line([xlim NaN I1x I1x],[I1y I1y NaN ylim],z,'color','m');
%      prompt  = {'X-coordinate at 1st corner:',...
%                 'Y-coordinate at 1st corner:'};
%      def     = {int2str(I1x),int2str(I1y)};
%      answer  = inputdlg(prompt,dlgTitle,lineNo,def);
%      answer  = str2num(char(answer{:}));
%      break
%     end
%    end
%    axes(ha)
%    
%    % Checks inputs:
%    if ~isempty(answer) && isfloat(answer) && (length(answer)==2) && ...
%                                                       all(isfinite(answer))
%     U1x = answer(1); U1y = answer(2);
%     secondcorner = true;
%    else
%     secondcorner = false;
%     warning('CVARGAS:ginput2:incorrectGeoreference',...
%             'Ignored incorrect georeference corners.')
%    end
%   
%    while secondcorner
%     % Selects second corner:
%     theans = ...
%          questdlg('Select the second corner (2 of 2):',dlgTitle,'OK','OK');
%     if ~strcmp(theans,'OK'), return, end
%     pause(secpause)
%    
%     [I2x,I2y] = ginput2(1,false,'none','UnZoom');
%     I2x       = round(I2x);
%     I2y       = round(I2y);
%     if ~ishandle(ha), return, end
%     if (ha==gca) && ~isempty(I2x) && ~isnan(I2x) && ...
%       (I2x~=I1x) && (I2y~=I1y)
%      axis(ha,[xlim ylim])
%      hgeo(2) = line([xlim NaN I2x I2x],[I2y I2y NaN ylim],z,'color','c');
%      prompt  = {'X-coordinate at 2nd corner:',...
%                 'Y-coordinate at 2nd corner:'};
%      def     = {int2str(I2x),int2str(I2y)};
%      answer  = inputdlg(prompt,dlgTitle,lineNo,def);
%      answer  = str2num(char(answer{:}));
%      break
%     end
%    end
%    axes(ha)
%    
%    % Checks inputs:
%    if secondcorner && ~isempty(answer) && isfloat(answer) && ...
%                          (length(answer)==2) && all(isfinite(answer))
%     U2x = answer(1); U2y = answer(2);
%    else
%     warning('CVARGAS:ginput2:incorrectGeoreference',...
%             'Ignored incorrect georeference corners.')
%    end
%   
%    % Deletes corner's lines:
%    if any(ishandle(hgeo))
%     delete(hgeo(ishandle(hgeo)))
%    end
%    
%    % Scale matrix:
%     SCALEMAT = [I1x I1y I2x I2y; U1x U1y U2x U2y];
%   else
%    % Continue
%   end
%  else
%   warning('CVARGAS:ginput2:noImageFound',...
%    'No image found in the current axes to georeference.')
%  end
%  
%  % OK, set the scaling then:
%  if ~isempty(SCALEMAT)
%   xdata = interp1(SCALEMAT(1,[1 3]),SCALEMAT(2,[1 3]),xdata,method,extrap);
%   ydata = interp1(SCALEMAT(1,[2 4]),SCALEMAT(2,[2 4]),ydata,method,extrap);
%   xlim2 = interp1(SCALEMAT(1,[1 3]),SCALEMAT(2,[1 3]),xlim ,method,extrap);
%   ylim2 = interp1(SCALEMAT(1,[2 4]),SCALEMAT(2,[2 4]),ylim ,method,extrap);
%   set(hi,'XData',xdata);
%   set(hi,'YData',ydata);  
%   set(ha,'XLim' ,sort(xlim2,'ascend'));
%   set(ha,'YLim' ,sort(ylim2,'ascend'));
%   % Reverses axis directions:
%   if diff(xlim)*diff(xlim2)<1
%    if strcmp(get(ha,'XDir'),'normal')
%     set(ha,'XDir','reverse')
%    else
%     set(ha,'XDir','normal')
%    end
%   end
%   if diff(ylim)*diff(ylim2)<1
%    if strcmp(get(ha,'YDir'),'normal')
%     set(ha,'YDir','reverse')
%    else
%     set(ha,'YDir','normal')
%    end
%   end
%  end
%  axis(ha,'normal')
%  
% end
% 
% % DRAWS?:
% if strcmpi(PlotOpt{1},'none')
%  yesdraw = false;
% else
%  yesdraw = true;
% end
% % Optional parameters:
% if yesdraw
%  hpoints  = [];
%  % Check for linestyle color:
%  yescolor     = true;
%  Nplotopt     = length(PlotOpt);
%  yeslinestyle = rem(Nplotopt,2);
%  if yeslinestyle % Given LineStyle
%   for k = 1:length(PlotOpt{1})
%    switch lower(PlotOpt{1}(k))
%     case 'y', yescolor = false; break
%     case 'm', yescolor = false; break
%     case 'c', yescolor = false; break
%     case 'r', yescolor = false; break
%     case 'g', yescolor = false; break
%     case 'b', yescolor = false; break
%     case 'w', yescolor = false; break
%     case 'k', yescolor = false; break
%     otherwise, % no color specified
%    end
%   end
%  end
%  if ~yescolor && (Nplotopt*yeslinestyle~=1)
%   for k = yeslinestyle+1:2:Nplotopt    % Given 'Color'
%    if strncmpi(PlotOpt{k},'co',2), yescolor = false; break, end
%   end
%  end
%  if yescolor
%   contnan  = 1;
%   colors   = get(gca,'ColorOrder');
%   ncolors  = size(colors,1);
%   color    = colors(1,:);
%  end
% end
% 
% 
% % -------------------------------------------------------------------------
% % MAIN
% % -------------------------------------------------------------------------
% 
% cont       = 0;
% alim.ha    = [];
% alim.la    = {};
% undoPtrFig = [];
% 
% while cont<N     % Principal loop
%  
%  % GINPUT:
%  try
%   [x,y,button] = ginput(1);
%  catch % Changed for compatibility.
%   % GINPUT error:
%   if YESERROR
%      error('CVARGAS:ginput2:executionError',lasterr)
%   else
%    warning('CVARGAS:ginput2:executionError',lasterr)
%    if nargout<2  % Fixed BUG 29 SEP, 2008
%     X = [X Y]; 
%    end
%    return
%   end
%  end
%  
%  % Axes clicked:
%  ha = gca;
%  
%  % Gets limits:
%  if ~any(alim.ha==ha)
%   alim.ha(end+1) = ha;
%   alim.la{end+1} = axis;
%  end
%  
%  % Sets zoom:
%  zlim = getappdata(ha,'zoom_zoomOrigAxesLimits');
%  if isempty(zlim) % Fixed BUG, SEP 2008
%   zoom reset
%   zlim = getappdata(ha,'zoom_zoomOrigAxesLimits');
%  end
%  
%  % Checks if DOUBLE clicks:
%  pause(secpause) % Gives time for response
%  if strcmp(get(gcf,'selectiontype'),'open')
%   button = DOUBLECLICK;
%  end
% 
%  % Checks if ENTER or FINISHES button:
%  if isempty(button) || ismember(button,FINISHES)
%   % Finishes selection:
%   if (N==1) && isempty(X) % New feature v3.1
%    BUTTON = button;
%   end
%   break
%  end
%  
%  % Checks if DELETION button:
%  if ismember(button,DELETES)
%   if ~isempty(X)
%    inan = isnan(X(end));
%    if yesdraw
%     if ~inan
%      % Deletes last drawn point:
%      if ~isempty(hpoints) && ishandle(hpoints(end)) % Fixed bug Aug 2009
%       delete(hpoints(end)), hpoints(end) = [];
%      end
%     elseif yescolor
%      % Change color as necessary:
%      contnan = contnan-1;
%      color   = colors(mod(contnan-1,ncolors)+1,:);
%     end
%    end
%    % Deletes the last selected point:
%    X(end)      = [];
%    Y(end)      = [];
%    BUTTON(end) = [];
%    % Checks if the last point was NaN:
%    if ~inan
%     cont = cont-1;
%    end
%   elseif N==1
%    % Finishes selection: New feature v3.1
%    BUTTON = button;
%    break
%   end
%   continue
%  end
%  
%  % Checks if ZOOM OUT button:
%  if ismember(button,ZOOMOUT)
%   undoPtrFig = gcf;
%   setptr(undoPtrFig,'glassminus')
%   zoom out
%   continue
%  end
%  
%  % Checks if ZOOM RESET button:
%  if ismember(button,ZOOMRESET)
%   zoom reset
%   continue
%  end
%  
%  % Checks if the mouse was inside an axes of the current figure;
%  lim     = axis;
%  outside = x<lim(1) || x>lim(2) || y<lim(3) || y>lim(4);
%  
%  % Checks if ZOOM IN with PAN:
%  if ismember(button,ZOOMIN) && ~outside
%   % Dragging rectangle:
%   undoPtrFig = gcf;
%   setptr(undoPtrFig,'closedhand')
%   rbbox
%   ydrag = get(gca,'CurrentPoint');
%   xdrag = ydrag(1,1)-x;
%   ydrag = ydrag(1,2)-y;
%   % Do the PANNING:
%   if  any(abs([xdrag ydrag])>eps*1000000)
%    % Only PAN (dragging):
%    lim(1:4) = lim(1:4) - [xdrag xdrag ydrag ydrag]; 
%    axis(lim)
%   else
%    % PAN (centers the last point) and ZOOM in:
%    setptr(undoPtrFig,'glassplus')
%    lim = [x+diff(lim(1:2))/2*[-1 1] y+diff(lim(3:4))/2*[-1 1]]; 
%    axis(lim)
%    zoom(2)
%   end
%   continue
%  end
%  
%  % Checks if SELECTS button:
%  if ismember(button,SELECTS)
%   
%   % Sets NaNs if outside the axes:
%   if outside 
%    if ~isnumeric(N) 
%     % Change color:
%     if yesdraw && yescolor 
%      contnan = contnan+1;
%      color   = colors(mod(contnan-1,ncolors)+1,:);
%     end
%     % Adds NaNs but the point counters do not take it into account:
%     X      = [X;      NaN];
%     Y      = [Y;      NaN];
%     BUTTON = [BUTTON; button];
%    else
%     % Ignores the point
%    end
%   else
%    % Draws the result:
%    if yesdraw
%     % Search for last point:
%     x0 = []; y0 = []; z0 = [];
%     inan = isnan([NaN; X; NaN]);
%     if ~inan(end-1)
%      inan        = find(inan);
%      nlastpoints = inan(end)-inan(end-1)-1;
%      npoints     = length(hpoints);
%      range       = npoints-nlastpoints+1:npoints;
%      hlastaxes   = get(hpoints(range),'Parent');
%      if iscell(hlastaxes), hlastaxes  = cell2mat(hlastaxes); end
%      [loc,loc] = ismember(ha,hlastaxes);
%      if loc
%       x0 = get(hpoints(range(loc)),'XData');
%       y0 = get(hpoints(range(loc)),'YData');
%       z0 = get(hpoints(range(loc)),'ZData');
%      end
%     end
%     holdon = ishold;
%     if ~holdon, hold on, end 
%      h = plot([x0 x],[y0 y],PlotOpt{:});
%      % Elevates the value:
%      z = get(ha,'ZLim'); z = z(2);
%      set(h,'Zdata',[z0 z])
%      % Sets the color:
%      if yescolor
%       set(h,'Color',color)
%      end
%      hpoints = [hpoints; h];
%     if ~holdon, hold off, end 
%     % Restores limits:
%     axis(lim)
%    end
%    
%    % Centers the selected point if ZOOM-IN: 29 SEP,2008
%    if all((lim~=zlim))
%     lim = [x+diff(lim(1:2))/2*[-1 1] y+diff(lim(3:4))/2*[-1 1]];
%     axis(lim)
%    end
%    
%    % Saves the result:
%    X      = [X;      x];
%    Y      = [Y;      y];
%    BUTTON = [BUTTON; button];
%    cont   = cont+1;
%   end
%   continue
%  end
% 
%  % Checks if any other button pressed inside the axes:
%  if ~outside
%   X      = [X;      x];
%   Y      = [Y;      y];
%   BUTTON = [BUTTON; button];
%   cont   = cont+1;
%  else
%   if N==1 % New feature v3.1
%    BUTTON = button;
%    break
%   end
%   % ignores the selection
%  end
%  
% end
% 
% % Returns pointer.
% if ~isempty(undoPtrFig) && ishandle(undoPtrFig)
%  setptr(undoPtrFig,'arrow')
% end
% 
% % Deletes drawn points if still exist:
% if yesdraw && any(ishandle(hpoints))
%  delete(hpoints(ishandle(hpoints)))
% end
% 
% % Returns original limits.
% if ~strcmp(UnZoom,'KeepZoom') && ~isempty(alim.ha)
%  alim.ha(~ishandle(alim.ha)) = [];
%  for k = 1:length(alim.ha)
%   temp = axis(alim.ha(k));
%   if ~all(temp(1:4)==alim.la{k}(1:4))
%    axis(alim.ha(k),alim.la{k})
%   end
%  end
% end
% 
% 
% % OUTPUTS CHECK-OUT
% % -------------------------------------------------------------------------
% 
% if nargout<2
%  X = [X Y]; 
% end
% 
% end
% % [EOF]   ginput2.m

 end
end

function [hdr, record] = edfRead(fname, varargin)
% Read European Data Format file into MATLAB
%
% [hdr, record] = edfRead(fname)
%         Reads data from ALL RECORDS of file fname ('*.edf'). Header
%         information is returned in structure hdr, and the signals
%         (waveforms) are returned in structure record, with waveforms
%         associated with the records returned as fields titled 'data' of
%         structure record.
%
% [...] = edfRead(fname, 'assignToVariables', assignToVariables)
%         Triggers writing of individual output variables, as defined by
%         field 'labels', into the caller workspace.
%
% FORMAT SPEC: Source: http://www.edfplus.info/specs/edf.html SEE ALSO:
% http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/eeg/edf_spec.htm
%
% The first 256 bytes of the header record specify the version number of
% this format, local patient and recording identification, time information
% about the recording, the number of data records and finally the number of
% signals (ns) in each data record. Then for each signal another 256 bytes
% follow in the header record, each specifying the type of signal (e.g.
% EEG, body temperature, etc.), amplitude calibration and the number of
% samples in each data record (from which the sampling frequency can be
% derived since the duration of a data record is also known). In this way,
% the format allows for different gains and sampling frequencies for each
% signal. The header record contains 256 + (ns * 256) bytes.
%
% Following the header record, each of the subsequent data records contains
% 'duration' seconds of 'ns' signals, with each signal being represented by
% the specified (in the header) number of samples. In order to reduce data
% size and adapt to commonly used software for acquisition, processing and
% graphical display of polygraphic signals, each sample value is
% represented as a 2-byte integer in 2's complement format. Figure 1 shows
% the detailed format of each data record.
%
% DATA SOURCE: Signals of various types (including the sample signal used
% below) are available from PHYSIONET: http://www.physionet.org/
%
%
% % EXAMPLE 1:
% % Read all waveforms/data associated with file 'ecgca998.edf':
%
% [header, recorddata] = edfRead('ecgca998.edf');
%
% % EXAMPLE 2:
% % Read records 3 and 5, associated with file 'ecgca998.edf':
%
% header = edfRead('ecgca998.edf','AssignToVariables',true);
% % Header file specifies data labels 'label_1'...'label_n'; these are
% % created as variables in the caller workspace.
%
% Coded 8/27/09 by Brett Shoelson, PhD
% brett.shoelson@mathworks.com
% Copyright 2009 - 2012 MathWorks, Inc.

% HEADER RECORD
% 8 ascii : version of this data format (0)
% 80 ascii : local patient identification
% 80 ascii : local recording identification
% 8 ascii : startdate of recording (dd.mm.yy)
% 8 ascii : starttime of recording (hh.mm.ss)
% 8 ascii : number of bytes in header record
% 44 ascii : reserved
% 8 ascii : number of data records (-1 if unknown)
% 8 ascii : duration of a data record, in seconds
% 4 ascii : number of signals (ns) in data record
% ns * 16 ascii : ns * label (e.g. EEG FpzCz or Body temp)
% ns * 80 ascii : ns * transducer type (e.g. AgAgCl electrode)
% ns * 8 ascii : ns * physical dimension (e.g. uV or degreeC)
% ns * 8 ascii : ns * physical minimum (e.g. -500 or 34)
% ns * 8 ascii : ns * physical maximum (e.g. 500 or 40)
% ns * 8 ascii : ns * digital minimum (e.g. -2048)
% ns * 8 ascii : ns * digital maximum (e.g. 2047)
% ns * 80 ascii : ns * prefiltering (e.g. HP:0.1Hz LP:75Hz)
% ns * 8 ascii : ns * nr of samples in each data record
% ns * 32 ascii : ns * reserved

% DATA RECORD
% nr of samples[1] * integer : first signal in the data record
% nr of samples[2] * integer : second signal
% ..
% ..
% nr of samples[ns] * integer : last signal

if nargin > 3
    error('EDFREAD: Too many input arguments.');
end

if ~nargin
    error('EDFREAD: Requires at least one input argument (filename to read).');
end

if nargin == 1
    assignToVariables = false;
end

[fid,msg] = fopen(fname,'r');
if fid == -1
    error(msg)
end

assignToVariables = false; %Default
for ii = 1:2:numel(varargin)
    switch lower(varargin{ii})
        case 'assigntovariables'
            assignToVariables = varargin{ii+1};
    end
end

% HEADER
hdr.ver        = str2double(char(fread(fid,8)'));
hdr.patientID  = fread(fid,80,'*char')';
hdr.recordID   = fread(fid,80,'*char')';
hdr.startdate  = fread(fid,8,'*char')';% (dd.mm.yy)
% hdr.startdate  = datestr(datenum(fread(fid,8,'*char')','dd.mm.yy'), 29); %'yyyy-mm-dd' (ISO 8601)
hdr.starttime  = fread(fid,8,'*char')';% (hh.mm.ss)
% hdr.starttime  = datestr(datenum(fread(fid,8,'*char')','hh.mm.ss'), 13); %'HH:MM:SS' (ISO 8601)
hdr.bytes      = str2double(fread(fid,8,'*char')');
reserved       = fread(fid,44);
hdr.records    = str2double(fread(fid,8,'*char')');
hdr.duration   = str2double(fread(fid,8,'*char')');
% Number of signals
hdr.ns    = str2double(fread(fid,4,'*char')');
for ii = 1:hdr.ns
    hdr.label{ii} = fread(fid,16,'*char')';
end
for ii = 1:hdr.ns
    hdr.transducer{ii} = fread(fid,80,'*char')';
end
% Physical dimension
for ii = 1:hdr.ns
    hdr.units{ii} = fread(fid,8,'*char')';
end
% Physical minimum
for ii = 1:hdr.ns
    hdr.physicalMin(ii) = str2double(fread(fid,8,'*char')');
end
% Physical maximum
for ii = 1:hdr.ns
    hdr.physicalMax(ii) = str2double(fread(fid,8,'*char')');
end
% Digital minimum
for ii = 1:hdr.ns
    hdr.digitalMin(ii) = str2double(fread(fid,8,'*char')');
end
% Digital maximum
for ii = 1:hdr.ns
    hdr.digitalMax(ii) = str2double(fread(fid,8,'*char')');
end
for ii = 1:hdr.ns
    hdr.prefilter{ii} = fread(fid,80,'*char')';
end
for ii = 1:hdr.ns
    hdr.samples(ii) = str2double(fread(fid,8,'*char')');
end
for ii = 1:hdr.ns
    reserved    = fread(fid,32,'*char')';
end
hdr.label = deblank(hdr.label);
hdr.units = deblank(hdr.units);


if nargout > 1 || assignToVariables
    % Scale data (linear scaling)
    scalefac = (hdr.physicalMax - hdr.physicalMin)./(hdr.digitalMax - hdr.digitalMin);
    dc = hdr.physicalMax - scalefac .* hdr.digitalMax;

    % RECORD DATA REQUESTED
    tmpdata = struct;
    for recnum = 1:hdr.records
        for ii = 1:hdr.ns
            % Use a cell array for DATA because number of samples may vary
            % from sample to sample
            tmpdata(recnum).data{ii} = fread(fid,hdr.samples(ii),'int16') * scalefac(ii) + dc(ii);
        end
    end
    record = zeros(hdr.ns, hdr.samples(1)*hdr.records);

    for ii = 1:numel(hdr.label)
        ctr = 1;
        for jj = 1:hdr.records
            try
                record(ii, ctr : ctr + hdr.samples - 1) = tmpdata(jj).data{ii};
            end
            ctr = ctr + hdr.samples;
        end
    end

    if assignToVariables
        for ii = 1:numel(hdr.label)
            try
                eval(['assignin(''caller'',''',hdr.label{ii},''',record(ii,:))'])
            end
        end
    end
end
fclose(fid);
end