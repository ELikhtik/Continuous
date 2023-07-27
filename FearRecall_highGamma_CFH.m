% using "pips" on Becky's data

%% 0 Initialization - Carolina-modified script from Becky, oct22 NON RMS DATA

clear all
close all
clc

Preprocessdir='C:\Users\ellab.GENECTR\Desktop\current_subjects\BECKYS_DATA\MERGED_FEMALES\RECALL\'; 
Postprocessdir = 'C:\Users\ellab.GENECTR\Desktop\current_subjects\BECKYS_DATA\MERGED_FEMALES\RECALL\processed_denoise60hz_filt70to150\';
dirname=Preprocessdir; 
Animaldir = dir([dirname,'TPH*']);
numCS=10; %Training CS = 5    Recall CS = 10    Context CS = 0
%load('IL_ephys_workspace.mat') %This will automatically load the data post-processing in the final structures for manipulation and figure creation. If you load this then skip to section 9!
A=dir([dirname,'TPH*']);

%% 1 Separating ephys data by trial and filtering
%Setting the directory and finding all animal files within the directory 

%This for loop will cut each electrode recording for each mouse into the trials and filter it
for animals=17
  
    clear LFPs
%     try
     filename=[dirname,Animaldir(animals).name]
      openNEV([filename,'\RECALL.nev'],'noread')
     openNSx([filename,'\RECALL.ns3']) %Use ns4 for Morgans data
    
    nCSCs=length(NS3.ElectrodesInfo);
    switch nCSCs
        case 2
            names={'CeA_L','BNST_L'};
       
    end
    
   NS3.Data=double(NS3.Data);
                 
    LFPs.ts=linspace(1,NS3.MetaTags.DataPoints,NS3.MetaTags.DataPoints);
    LFPs.filename=Animaldir(animals).name;
                     
        LFPs.dn60preRMS=struct('CeA_L',[],'BNST_L',[]);  %data has 60hz noise
        
        for a=1:length(names)
        LFPs.dn60preRMS.(names{a})=removeLineNoise_SpectrumEstimation(NS3.Data(a,:),2000,['NH=5','LF=60']);
        end          
   
    trialstart=double(NEV.Data.SerialDigitalIO.TimeStampSec(NEV.Data.SerialDigitalIO.UnparsedData==65503)); %trialstart points come from NEV file
    trialstop=NEV.Data.SerialDigitalIO.TimeStampSec(NEV.Data.SerialDigitalIO.UnparsedData==65535);
    LFPs.OnTimes=trialstart-1; 
    LFPs.OffTimes=trialstart + 31;  %LFPs.OffTimes=trialstop;
    
   for a=1:length(LFPs.OnTimes)
        start = LFPs.OnTimes(a)*2000;
        precsStart = (LFPs.OnTimes(a)-30)*2000;
        stop = LFPs.OffTimes(a)*2000;
        precsStop = (LFPs.OffTimes(a)-30)*2000;
        LFPs.CSts{a}= LFPs.ts(1,start:stop-1);
        LFPs.pCSts{a}=LFPs.ts(1,precsStart:precsStop-1);
   end  
    
  for b=1:length(names)
        for c=1:length(LFPs.CSts)
            ind=cell2mat(LFPs.CSts(c));
            ind(1);
            ind(length(ind));
            [LFPs.Filt.(names{b}).CS{c},foo,LFPs.Pwr.(names{b}).CS{c},LFPs.Phase.(names{b}).CS{c}]=SinoFilt_Joe((LFPs.dn60preRMS.(names{b})(ind(1):ind(length(ind)))),2000,8,70,150);
            end
        end   
  
     for b=1:length(names)
        for c=1:length(LFPs.pCSts)
            ind=cell2mat(LFPs.pCSts(c));
            ind(1);
            ind(length(ind));
            [LFPs.Filt.(names{b}).pCS{c},foo,LFPs.Pwr.(names{b}).pCS{c},LFPs.Phase.(names{b}).pCS{c}]=SinoFilt_Joe((LFPs.dn60preRMS.(names{b})(ind(1):ind(length(ind)))),2000,8,70,150);
              end
    end 
        
 for b=1:length(names)
        for c=1:length(LFPs.CSts)
            ind=cell2mat(LFPs.CSts(c));
            ind(1);
            ind(length(ind));
            LFPs.unFilt.(names{b}).CS{c}=LFPs.dn60preRMS.(names{b})(ind(1):ind(length(ind)));
            end
        end   
  
     for b=1:length(names)
        for c=1:length(LFPs.pCSts)
            ind=cell2mat(LFPs.pCSts(c));
            ind(1);
            ind(length(ind));
            LFPs.unFilt.(names{b}).pCS{c}=LFPs.dn60preRMS.(names{b})(ind(1):ind(length(ind)));
              end
    end


         save([Postprocessdir,LFPs.filename],'LFPs','-v7.3');
     clear LFPs
end

%% 2 Sorting animals into experimental(TPH+ = TPHPOS) and control(TPH- = TPHNEG)
%Here we change the directory into the newly created processed folder
%containing all the cut and filtered data
dirname=Postprocessdir;

LFPsall.TPHNEG=struct('TPH645_RECALL',[],'TPH646_RECALL',[],'TPH647_RECALL',[],'TPH662_RECALL',[],'TPH891_RECALL',[],'TPH906_RECALL',[],'TPH1006_RECALL',[],'TPH1025_RECALL',[],'TPH1398_RECALL',[],'TPH1406_RECALL',[],'TPH1474_RECALL',[],'TPH1480_RECALL',[],'TPH1484_RECALL',[],'TPH1486_RECALL',[],'TPH1509_RECALL',[]);

LFPsall.TPHPOS=struct('TPH394_RECALL',[],'TPH655_RECALL',[],'TPH657_RECALL',[],'TPH691_RECALL',[],'TPH905_RECALL',[],'TPH907_RECALL',[],'TPH1008_RECALL',[],'TPH1026_RECALL',[],'TPH1389_RECALL',[],'TPH1390_RECALL',[],'TPH1397_RECALL',[],'TPH1399_RECALL',[],'TPH1477_RECALL',[],'TPH1487_RECALL',[]);


%This for loop splits the animals into experimental or control groups and 
%places them into the structure we will be pulling from for the rest of the code.
Animaldir=dir([dirname, 'TPH*']); 
for i=1:length(Animaldir)
    A=dir([dirname,Animaldir(i).name]);
    filename=[dirname,Animaldir(i).name];
    load(filename)
    anim = i;
    switch anim
        case 1
        LFPsall.TPHNEG.(['TPH645_RECALL'])=LFPs;
        case 2
        LFPsall.TPHNEG.(['TPH646_RECALL'])=LFPs;
        case 3
        LFPsall.TPHNEG.(['TPH647_RECALL'])=LFPs;
        case 4
        LFPsall.TPHNEG.(['TPH662_RECALL'])=LFPs;
        case 5
        LFPsall.TPHNEG.(['TPH891_RECALL'])=LFPs;
        case 6
        LFPsall.TPHNEG.(['TPH906_RECALL'])=LFPs;
        case 7
        LFPsall.TPHNEG.(['TPH1006_RECALL'])=LFPs;
        case 8
        LFPsall.TPHNEG.(['TPH1025_RECALL'])=LFPs;
        case 9
        LFPsall.TPHNEG.(['TPH1398_RECALL'])=LFPs;
        case 10
        LFPsall.TPHNEG.(['TPH1406_RECALL'])=LFPs;  
        case 11
        LFPsall.TPHNEG.(['TPH1509_RECALL'])=LFPs;
        case 12
        LFPsall.TPHPOS.(['TPH394_RECALL'])=LFPs;
        case 13
        LFPsall.TPHPOS.(['TPH655_RECALL'])=LFPs;
        case 14
        LFPsall.TPHPOS.(['TPH657_RECALL'])=LFPs;
        case 15
        LFPsall.TPHPOS.(['TPH691_RECALL'])=LFPs;
        case 16
        LFPsall.TPHPOS.(['TPH905_RECALL'])=LFPs;
        case 17
        LFPsall.TPHPOS.(['TPH907_RECALL'])=LFPs;
        case 18
        LFPsall.TPHPOS.(['TPH1008_RECALL'])=LFPs;
        case 19
        LFPsall.TPHPOS.(['TPH1026_RECALL'])=LFPs;
        case 20
        LFPsall.TPHPOS.(['TPH1389_RECALL'])=LFPs;
        case 22
        LFPsall.TPHPOS.(['TPH1390_RECALL'])=LFPs;
        case 23
        LFPsall.TPHPOS.(['TPH1397_RECALL'])=LFPs;
        case 24
        LFPsall.TPHPOS.(['TPH1399_RECALL'])=LFPs;
                      
    end
        clear LFPs
end

%% 3 Power Calculation with mtcsg for TPHNEG & TPHPOS, graphing frequency power, area under curve calculation

%TPHNEG TonePwr Calculation CS & pCS
names={'CeA_L','BNST_L'};

TPHNEG={'TPH645_RECALL','TPH646_RECALL','TPH647_RECALL','TPH662_RECALL','TPH891_RECALL','TPH906_RECALL','TPH1006_RECALL','TPH1025_RECALL','TPH1398_RECALL','TPH1406_RECALL','TPH1509_RECALL',};

TPHPOS={'TPH394_RECALL','TPH655_RECALL','TPH657_RECALL','TPH691_RECALL','TPH905_RECALL','TPH907_RECALL','TPH1008_RECALL','TPH1026_RECALL','TPH1389_RECALL','TPH1390_RECALL','TPH1397_RECALL','TPH1399_RECALL'};

CSname={'CS1','CS2','CS3','CS4','CS5','CS6','CS7','CS8','CS9','CS10'};
pCSname={'pCS1','pCS2','pCS3','pCS4','pCS5','pCS6','pCS7','pCS8','pCS9','pCS10'};

Pipname={'pip1','pip2','pip3','pip4','pip5','pip6','pip7','pip8','pip9','pip10','pip11','pip12','pip13','pip14','pip15','pip16','pip17','pip18','pip19','pip20','pip21','pip22','pip23','pip24','pip25','pip26','pip27','pip28','pip29','pip30'};
Prepipname={'prepip1','prepip2','prepip3','prepip4','prepip5','prepip6','prepip7','prepip8','prepip9','prepip10','prepip11','prepip12','prepip13','prepip14','prepip15','prepip16','prepip17','prepip18','prepip19','prepip20','prepip21','prepip22','prepip23','prepip24','prepip25','prepip26','prepip27','prepip28','prepip29','prepip30'};

for a=1:length(TPHNEG)
    for b=1:length(names)   
       LFPsall.TPHNEG.(TPHNEG{a}).TonePwr.(names{b}).CSs=struct('CS1',[],'CS2',[],'CS3',[],'CS4',[],'CS5',[],'CS6',[],'CS7',[],'CS8',[],'CS9',[],'CS10',[]);
       LFPsall.TPHNEG.(TPHNEG{a}).TonePwr.(names{b}).pCSs=struct('pCS1',[],'pCS2',[],'pCS3',[],'pCS4',[],'pCS5',[],'pCS6',[],'pCS7',[],'pCS8',[],'pCS9',[],'pCS10',[]);
       LFPsall.TPHNEG.(TPHNEG{a}).TonePwr.(names{b}).pips.Pips_freqtime=struct('CS1',[],'CS2',[],'CS3',[],'CS4',[],'CS5',[],'CS6',[],'CS7',[],'CS8',[],'CS9',[],'CS10',[],'all_CSs',[], 'all_CSs_mean',[], 'CS1_5',[], 'CS1_5_mean',[], 'CS6_10',[], 'CS6_10_mean',[]);
       LFPsall.TPHNEG.(TPHNEG{a}).TonePwr.(names{b}).pips.Prepips_freqtime=struct('pCS1',[],'pCS2',[],'pCS3',[],'pCS4',[],'pCS5',[],'pCS6',[],'pCS7',[],'pCS8',[],'pCS9',[],'pCS10',[], 'all_CSs',[], 'all_CSs_mean',[], 'pCS1_5',[], 'pCS1_5_mean',[], 'pCS6_10',[], 'pCS6_10_mean',[]);
  
    end
end

for a=1:length(TPHNEG)
    for b=1:length(names)
        for c=1:numCS            
        LFPsall.TPHNEG.(TPHNEG{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c})=struct('pip1',[],'pip2',[],'pip3',[],'pip4',[],'pip5',[],'pip6',[],'pip7',[],'pip8',[],'pip9',[],'pip10',[],'pip11',[],'pip12',[],'pip13',[],'pip14',[],'pip15',[],'pip16',[],'pip17',[],'pip18',[],'pip19',[],'pip20',[],'pip21',[],'pip22',[],'pip23',[],'pip24',[],'pip25',[],'pip26',[],'pip27',[],'pip28',[],'pip29',[],'pip30',[],'altogether',[],'altogether_mean',[],'pip1_mean',[],'pip2_mean',[],'pip3_mean',[],'pip4_mean',[],'pip5_mean',[],'pip6_mean',[],'pip7_mean',[],'pip8_mean',[],'pip9_mean',[],'pip10_mean',[],'pip11_mean',[],'pip12_mean',[],'pip13_mean',[],'pip14_mean',[],'pip15_mean',[],'pip16_mean',[],'pip17_mean',[],'pip18_mean',[],'pip19_mean',[],'pip20_mean',[],'pip21_mean',[],'pip22_mean',[],'pip23_mean',[],'pip24_mean',[],'pip25_mean',[],'pip26_mean',[],'pip27_mean',[],'pip28_mean',[],'pip29_mean',[],'pip30_mean',[],'meaneachpip_altogether',[]);
        LFPsall.TPHNEG.(TPHNEG{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c})=struct('prepip1',[],'prepip2',[],'prepip3',[],'prepip4',[],'prepip5',[],'prepip6',[],'prepip7',[],'prepip8',[],'prepip9',[],'prepip10',[],'prepip11',[],'prepip12',[],'prepip13',[],'prepip14',[],'prepip15',[],'prepip16',[],'prepip17',[],'prepip18',[],'prepip19',[],'prepip20',[],'prepip21',[],'prepip22',[],'prepip23',[],'prepip24',[],'prepip25',[],'prepip26',[],'prepip27',[],'prepip28',[],'prepip29',[],'prepip30',[],'altogether',[], 'altogether_mean',[],'prepip1_mean',[],'prepip2_mean',[],'prepip3_mean',[],'prepip4_mean',[],'prepip5_mean',[],'prepip6_mean',[],'prepip7_mean',[],'prepip8_mean',[],'prepip9_mean',[],'prepip10_mean',[],'prepip11_mean',[],'prepip12_mean',[],'prepip13_mean',[],'prepip14_mean',[],'prepip15_mean',[],'prepip16_mean',[],'prepip17_mean',[],'prepip18_mean',[],'prepip19_mean',[],'prepip20_mean',[],'prepip21_mean',[],'prepip22_mean',[],'prepip23_mean',[],'prepip24_mean',[],'prepip25_mean',[],'prepip26_mean',[],'prepip27_mean',[],'prepip28_mean',[],'prepip29_mean',[],'prepip30_mean',[],'meaneachpip_altogether',[]);
        end
    end
end

%convolution for tones
for a=1:length(TPHNEG)
    for b=1:length(names)
        for c=1:numCS
            [LFPsall.TPHNEG.(TPHNEG{a}).TonePwr.(names{b}).CSs.(CSname{c}),Freq_pw,Time_pw]=mtcsg(LFPsall.TPHNEG.(TPHNEG{a}).Filt.(names{b}).CS{1, c},2048,2000,500,480,2); 
            [LFPsall.TPHNEG.(TPHNEG{a}).TonePwr.(names{b}).pCSs.(pCSname{c}),Freq_pw,Time_pw]=mtcsg(LFPsall.TPHNEG.(TPHNEG{a}).Filt.(names{b}).pCS{1, c},2048,2000,500,480,2); 
        end
    end
end

%tphpos 

for a=1:length(TPHPOS)
    for b=1:length(names)   
       LFPsall.TPHPOS.(TPHPOS{a}).TonePwr.(names{b}).CSs=struct('CS1',[],'CS2',[],'CS3',[],'CS4',[],'CS5',[],'CS6',[],'CS7',[],'CS8',[],'CS9',[],'CS10',[]);
       LFPsall.TPHPOS.(TPHPOS{a}).TonePwr.(names{b}).pCSs=struct('pCS1',[],'pCS2',[],'pCS3',[],'pCS4',[],'pCS5',[],'pCS6',[],'pCS7',[],'pCS8',[],'pCS9',[],'pCS10',[]);
       LFPsall.TPHPOS.(TPHPOS{a}).TonePwr.(names{b}).pips.Pips_freqtime=struct('CS1',[],'CS2',[],'CS3',[],'CS4',[],'CS5',[],'CS6',[],'CS7',[],'CS8',[],'CS9',[],'CS10',[],'all_CSs',[], 'all_CSs_mean',[], 'CS1_5',[], 'CS1_5_mean',[], 'CS6_10',[], 'CS6_10_mean',[]);
       LFPsall.TPHPOS.(TPHPOS{a}).TonePwr.(names{b}).pips.Prepips_freqtime=struct('pCS1',[],'pCS2',[],'pCS3',[],'pCS4',[],'pCS5',[],'pCS6',[],'pCS7',[],'pCS8',[],'pCS9',[],'pCS10',[], 'all_CSs',[], 'all_CSs_mean',[], 'pCS1_5',[], 'pCS1_5_mean',[], 'pCS6_10',[], 'pCS6_10_mean',[]);
  
    end
end

for a=1:length(TPHPOS)
    for b=1:length(names)
        for c=1:numCS            
        LFPsall.TPHPOS.(TPHPOS{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c})=struct('pip1',[],'pip2',[],'pip3',[],'pip4',[],'pip5',[],'pip6',[],'pip7',[],'pip8',[],'pip9',[],'pip10',[],'pip11',[],'pip12',[],'pip13',[],'pip14',[],'pip15',[],'pip16',[],'pip17',[],'pip18',[],'pip19',[],'pip20',[],'pip21',[],'pip22',[],'pip23',[],'pip24',[],'pip25',[],'pip26',[],'pip27',[],'pip28',[],'pip29',[],'pip30',[],'altogether',[],'altogether_mean',[],'pip1_mean',[],'pip2_mean',[],'pip3_mean',[],'pip4_mean',[],'pip5_mean',[],'pip6_mean',[],'pip7_mean',[],'pip8_mean',[],'pip9_mean',[],'pip10_mean',[],'pip11_mean',[],'pip12_mean',[],'pip13_mean',[],'pip14_mean',[],'pip15_mean',[],'pip16_mean',[],'pip17_mean',[],'pip18_mean',[],'pip19_mean',[],'pip20_mean',[],'pip21_mean',[],'pip22_mean',[],'pip23_mean',[],'pip24_mean',[],'pip25_mean',[],'pip26_mean',[],'pip27_mean',[],'pip28_mean',[],'pip29_mean',[],'pip30_mean',[],'meaneachpip_altogether',[]);
        LFPsall.TPHPOS.(TPHPOS{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c})=struct('prepip1',[],'prepip2',[],'prepip3',[],'prepip4',[],'prepip5',[],'prepip6',[],'prepip7',[],'prepip8',[],'prepip9',[],'prepip10',[],'prepip11',[],'prepip12',[],'prepip13',[],'prepip14',[],'prepip15',[],'prepip16',[],'prepip17',[],'prepip18',[],'prepip19',[],'prepip20',[],'prepip21',[],'prepip22',[],'prepip23',[],'prepip24',[],'prepip25',[],'prepip26',[],'prepip27',[],'prepip28',[],'prepip29',[],'prepip30',[],'altogether',[], 'altogether_mean',[],'prepip1_mean',[],'prepip2_mean',[],'prepip3_mean',[],'prepip4_mean',[],'prepip5_mean',[],'prepip6_mean',[],'prepip7_mean',[],'prepip8_mean',[],'prepip9_mean',[],'prepip10_mean',[],'prepip11_mean',[],'prepip12_mean',[],'prepip13_mean',[],'prepip14_mean',[],'prepip15_mean',[],'prepip16_mean',[],'prepip17_mean',[],'prepip18_mean',[],'prepip19_mean',[],'prepip20_mean',[],'prepip21_mean',[],'prepip22_mean',[],'prepip23_mean',[],'prepip24_mean',[],'prepip25_mean',[],'prepip26_mean',[],'prepip27_mean',[],'prepip28_mean',[],'prepip29_mean',[],'prepip30_mean',[],'meaneachpip_altogether',[]);
        end
    end
end

%convolution for tones
for a=1:length(TPHPOS)
    for b=1:length(names)
        for c=1:numCS
            [LFPsall.TPHPOS.(TPHPOS{a}).TonePwr.(names{b}).CSs.(CSname{c}),Freq_pw,Time_pw]=mtcsg(LFPsall.TPHPOS.(TPHPOS{a}).Filt.(names{b}).CS{1, c},2048,2000,500,480,2);
            [LFPsall.TPHPOS.(TPHPOS{a}).TonePwr.(names{b}).pCSs.(pCSname{c}),Freq_pw,Time_pw]=mtcsg(LFPsall.TPHPOS.(TPHPOS{a}).Filt.(names{b}).pCS{1, c},2048,2000,500,480,2);
        end
    end
end

%% 4 Cut for pips
%we will now cut in a way that leaves each pip with 100bins
num_pips=30;
num_prepips=30;
pips_bins=zeros(num_pips,100);
prepips_bins=zeros(num_prepips,100); 

%preCS_bins=zeros(numCS,500);
%preCS_bins_spctg=zeros(numCS,100);
    
C=0;
    m=0;   
for n=1:num_pips  
    pips_bins(n,:)=88 + C.*m:1:88+99 + C.*m; 
    C=100;
    m=m+1;
end

C=0;
    m=0;   
for n=1:num_prepips  
    prepips_bins(n,:)=88 + C.*m:1:88+99 + C.*m; 
    C=100;
    m=m+1;
end


%now we isolate the power for the bins corresponding to each pip
 for a=1:length(TPHNEG)
     for b=1:length(names)
        for c=1:numCS
            for d=1:num_pips     
                 LFPsall.TPHNEG.(TPHNEG{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d})=(LFPsall.TPHNEG.(TPHNEG{a}).TonePwr.(names{b}).CSs.(CSname{c})(:,pips_bins(d,1):pips_bins(d,100))); %100 bins is 1000ms
             end
        end
     end
 end  
 
 %now we isolate the power for the bins corresponding to each prepip
  for a=1:length(TPHNEG)
     for b=1:length(names)
        for c=1:numCS
            for d=1:num_prepips     
                 LFPsall.TPHNEG.(TPHNEG{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d})=(LFPsall.TPHNEG.(TPHNEG{a}).TonePwr.(names{b}).pCSs.(pCSname{c})(:,prepips_bins(d,1):prepips_bins(d,100)));
             end
        end
     end
 end 

 for a=1:length(TPHPOS)
     for b=1:length(names)
        for c=1:numCS
            for d=1:num_pips     
                 LFPsall.TPHPOS.(TPHPOS{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d})=(LFPsall.TPHPOS.(TPHPOS{a}).TonePwr.(names{b}).CSs.(CSname{c})(:,pips_bins(d,1):pips_bins(d,100))); %100 bins is 1000ms
             end
        end
     end
 end  
 
 %now we isolate the power for the bins corresponding to each prepip
  for a=1:length(TPHPOS)
     for b=1:length(names)
        for c=1:numCS
            for d=1:num_prepips     
                 LFPsall.TPHPOS.(TPHPOS{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d})=(LFPsall.TPHPOS.(TPHPOS{a}).TonePwr.(names{b}).pCSs.(pCSname{c})(:,prepips_bins(d,1):prepips_bins(d,100)));
             end
        end
     end
  end 

  %% select good animals (placements only)

TPHNEG_CeA_L={'TPH645_RECALL','TPH646_RECALL','TPH647_RECALL','TPH891_RECALL','TPH1006_RECALL','TPH1406_RECALL','TPH1509_RECALL'};
TPHNEG_BNST_L={'TPH645_RECALL','TPH646_RECALL','TPH647_RECALL','TPH662_RECALL','TPH891_RECALL','TPH906_RECALL','TPH1006_RECALL','TPH1025_RECALL','TPH1398_RECALL','TPH1406_RECALL'};

TPHPOS_CeA_L={'TPH655_RECALL','TPH691_RECALL','TPH905_RECALL','TPH907_RECALL','TPH1008_RECALL','TPH1026_RECALL','TPH1390_RECALL'};
TPHPOS_BNST_L={'TPH394_RECALL','TPH655_RECALL','TPH657_RECALL','TPH905_RECALL','TPH907_RECALL','TPH1008_RECALL','TPH1026_RECALL','TPH1389_RECALL','TPH1390_RECALL','TPH1397_RECALL','TPH1399_RECALL'};

%% 6 Power calculations

%TPHNEG
%BNST

for a=1:length(TPHNEG_BNST_L)
    for b=2%:length(names)
        for c=1:numCS            
            ooooo=LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).altogether=cat(2,ooooo.pip1,ooooo.pip2,ooooo.pip3,ooooo.pip4,ooooo.pip5,ooooo.pip6,ooooo.pip7,ooooo.pip8,ooooo.pip9,ooooo.pip10,ooooo.pip11,ooooo.pip12,ooooo.pip13,ooooo.pip14,ooooo.pip15,ooooo.pip16,ooooo.pip17,ooooo.pip18,ooooo.pip19,ooooo.pip20,ooooo.pip21,ooooo.pip22,ooooo.pip23,ooooo.pip24,ooooo.pip25,ooooo.pip26,ooooo.pip27,ooooo.pip28,ooooo.pip29,ooooo.pip30);
        
            ooooo=LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).altogether=cat(2,ooooo.prepip1,ooooo.prepip2,ooooo.prepip3,ooooo.prepip4,ooooo.prepip5,ooooo.prepip6,ooooo.prepip7,ooooo.prepip8,ooooo.prepip9,ooooo.prepip10,ooooo.prepip11,ooooo.prepip12,ooooo.prepip13,ooooo.prepip14,ooooo.prepip15,ooooo.prepip16,ooooo.prepip17,ooooo.prepip18,ooooo.prepip19,ooooo.prepip20,ooooo.prepip21,ooooo.prepip22,ooooo.prepip23,ooooo.prepip24,ooooo.prepip25,ooooo.prepip26,ooooo.prepip27,ooooo.prepip28,ooooo.prepip29,ooooo.prepip30);
        end
    end
end

for a=1:length(TPHNEG_BNST_L)
    for b=2%:length(names) 
    ooooo=LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime;
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether,ooooo.CS5.altogether,ooooo.CS6.altogether,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs_mean=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs,2,'omitnan');
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4_mean=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4,2,'omitnan');
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10=cat(2,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10_mean=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10,2,'omitnan');
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2_mean=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10.altogether,2);
    end
end

for a=1:length(TPHNEG_BNST_L)
    for b=2%:length(names) 
    ooooo=LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime;
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs=cat(2,ooooo.pCS1.altogether,ooooo.pCS2.altogether,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs_mean=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs,2,'omitnan');
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10=cat(2,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10_mean=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10,2,'omitnan');     
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10.altogether,2);
    end
end

pip_mean_name={'pip1_mean','pip2_mean','pip3_mean','pip4_mean','pip5_mean','pip6_mean','pip7_mean','pip8_mean','pip9_mean','pip10_mean','pip11_mean','pip12_mean','pip13_mean','pip14_mean','pip15_mean','pip16_mean','pip17_mean','pip18_mean','pip19_mean','pip20_mean','pip21_mean','pip22_mean','pip23_mean','pip24_mean','pip25_mean','pip26_mean','pip27_mean','pip28_mean','pip29_mean','pip30_mean'};
prepip_mean_name={'prepip1_mean','prepip2_mean','prepip3_mean','prepip4_mean','prepip5_mean','prepip6_mean','prepip7_mean','prepip8_mean','prepip9_mean','prepip10_mean','prepip11_mean','prepip12_mean','prepip13_mean','prepip14_mean','prepip15_mean','prepip16_mean','prepip17_mean','prepip18_mean','prepip19_mean','prepip20_mean','prepip21_mean','prepip22_mean','prepip23_mean','prepip24_mean','prepip25_mean','prepip26_mean','prepip27_mean','prepip28_mean','prepip29_mean','prepip30_mean'};

for a=1:length(TPHNEG_BNST_L)
    for b=2%:length(names) 
        for c=1:length(CSname) 
            for d=1:length(Pipname) 
        LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(pip_mean_name{d})=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d}),2,'omitnan');
         end
    end
    end
end


for a=1:length(TPHNEG_BNST_L)
    for b=2%:length(names) 
        for c=1:length(pCSname) 
            for d=1:length(Prepipname) 
        LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(prepip_mean_name{d})=mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d}),2,'omitnan');
         end
    end
    end
end

for a=1:length(TPHNEG_BNST_L)
    for b=2%:length(names) 
        for c=1:length(CSname) 
            ooooo=LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.pip1_mean,ooooo.pip2_mean,ooooo.pip3_mean,ooooo.pip4_mean,ooooo.pip5_mean,ooooo.pip6_mean,ooooo.pip7_mean,ooooo.pip8_mean,ooooo.pip9_mean,ooooo.pip10_mean,ooooo.pip11_mean,ooooo.pip12_mean,ooooo.pip13_mean,ooooo.pip14_mean,ooooo.pip15_mean,ooooo.pip16_mean,ooooo.pip17_mean,ooooo.pip18_mean,ooooo.pip19_mean,ooooo.pip20_mean,ooooo.pip21_mean,ooooo.pip22_mean,ooooo.pip23_mean,ooooo.pip24_mean,ooooo.pip25_mean,ooooo.pip26_mean,ooooo.pip27_mean,ooooo.pip28_mean,ooooo.pip29_mean,ooooo.pip30_mean);  
        end
    end
end

for a=1:length(TPHNEG_BNST_L)
    for b=2%:length(names) 
        for c=1:length(pCSname) 
              ooooo=LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.prepip1_mean,ooooo.prepip2_mean,ooooo.prepip3_mean,ooooo.prepip4_mean,ooooo.prepip5_mean,ooooo.prepip6_mean,ooooo.prepip7_mean,ooooo.prepip8_mean,ooooo.prepip9_mean,ooooo.prepip10_mean,ooooo.prepip11_mean,ooooo.prepip12_mean,ooooo.prepip13_mean,ooooo.prepip14_mean,ooooo.prepip15_mean,ooooo.prepip16_mean,ooooo.prepip17_mean,ooooo.prepip18_mean,ooooo.prepip19_mean,ooooo.prepip20_mean,ooooo.prepip21_mean,ooooo.prepip22_mean,ooooo.prepip23_mean,ooooo.prepip24_mean,ooooo.prepip25_mean,ooooo.prepip26_mean,ooooo.prepip27_mean,ooooo.prepip28_mean,ooooo.prepip29_mean,ooooo.prepip30_mean);  
        end
    end
end

%TPHPOS
%BNST

for a=1:length(TPHPOS_BNST_L)
    for b=2%:length(names)
        for c=1:numCS            
            ooooo=LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).altogether=cat(2,ooooo.pip1,ooooo.pip2,ooooo.pip3,ooooo.pip4,ooooo.pip5,ooooo.pip6,ooooo.pip7,ooooo.pip8,ooooo.pip9,ooooo.pip10,ooooo.pip11,ooooo.pip12,ooooo.pip13,ooooo.pip14,ooooo.pip15,ooooo.pip16,ooooo.pip17,ooooo.pip18,ooooo.pip19,ooooo.pip20,ooooo.pip21,ooooo.pip22,ooooo.pip23,ooooo.pip24,ooooo.pip25,ooooo.pip26,ooooo.pip27,ooooo.pip28,ooooo.pip29,ooooo.pip30);
        
            ooooo=LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).altogether=cat(2,ooooo.prepip1,ooooo.prepip2,ooooo.prepip3,ooooo.prepip4,ooooo.prepip5,ooooo.prepip6,ooooo.prepip7,ooooo.prepip8,ooooo.prepip9,ooooo.prepip10,ooooo.prepip11,ooooo.prepip12,ooooo.prepip13,ooooo.prepip14,ooooo.prepip15,ooooo.prepip16,ooooo.prepip17,ooooo.prepip18,ooooo.prepip19,ooooo.prepip20,ooooo.prepip21,ooooo.prepip22,ooooo.prepip23,ooooo.prepip24,ooooo.prepip25,ooooo.prepip26,ooooo.prepip27,ooooo.prepip28,ooooo.prepip29,ooooo.prepip30);
        end
    end
end

for a=1:length(TPHPOS_BNST_L)
    for b=2%:length(names) 
    ooooo=LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime;
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether,ooooo.CS5.altogether,ooooo.CS6.altogether,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs_mean=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs,2,'omitnan');
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4_mean=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4,2,'omitnan');
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10=cat(2,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10_mean=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10,2,'omitnan');
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2_mean=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10.altogether,2);
    end
end


for a=1:length(TPHPOS_BNST_L)
    for b=2%:length(names) 
    ooooo=LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime;
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs=cat(2,ooooo.pCS1.altogether,ooooo.pCS2.altogether,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs_mean=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs,2,'omitnan');
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10=cat(2,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10_mean=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10,2,'omitnan');     
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10.altogether,2);
    end
end

pip_mean_name={'pip1_mean','pip2_mean','pip3_mean','pip4_mean','pip5_mean','pip6_mean','pip7_mean','pip8_mean','pip9_mean','pip10_mean','pip11_mean','pip12_mean','pip13_mean','pip14_mean','pip15_mean','pip16_mean','pip17_mean','pip18_mean','pip19_mean','pip20_mean','pip21_mean','pip22_mean','pip23_mean','pip24_mean','pip25_mean','pip26_mean','pip27_mean','pip28_mean','pip29_mean','pip30_mean'};
prepip_mean_name={'prepip1_mean','prepip2_mean','prepip3_mean','prepip4_mean','prepip5_mean','prepip6_mean','prepip7_mean','prepip8_mean','prepip9_mean','prepip10_mean','prepip11_mean','prepip12_mean','prepip13_mean','prepip14_mean','prepip15_mean','prepip16_mean','prepip17_mean','prepip18_mean','prepip19_mean','prepip20_mean','prepip21_mean','prepip22_mean','prepip23_mean','prepip24_mean','prepip25_mean','prepip26_mean','prepip27_mean','prepip28_mean','prepip29_mean','prepip30_mean'};

for a=1:length(TPHPOS_BNST_L)
    for b=2%:length(names) 
        for c=1:length(CSname) 
            for d=1:length(Pipname) 
        LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(pip_mean_name{d})=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d}),2,'omitnan');
         end
    end
    end
end


for a=1:length(TPHPOS_BNST_L)
    for b=2%:length(names) 
        for c=1:length(pCSname) 
            for d=1:length(Prepipname) 
        LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(prepip_mean_name{d})=mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d}),2,'omitnan');
         end
    end
    end
end

for a=1:length(TPHPOS_BNST_L)
    for b=2%:length(names) 
        for c=1:length(CSname) 
            ooooo=LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.pip1_mean,ooooo.pip2_mean,ooooo.pip3_mean,ooooo.pip4_mean,ooooo.pip5_mean,ooooo.pip6_mean,ooooo.pip7_mean,ooooo.pip8_mean,ooooo.pip9_mean,ooooo.pip10_mean,ooooo.pip11_mean,ooooo.pip12_mean,ooooo.pip13_mean,ooooo.pip14_mean,ooooo.pip15_mean,ooooo.pip16_mean,ooooo.pip17_mean,ooooo.pip18_mean,ooooo.pip19_mean,ooooo.pip20_mean,ooooo.pip21_mean,ooooo.pip22_mean,ooooo.pip23_mean,ooooo.pip24_mean,ooooo.pip25_mean,ooooo.pip26_mean,ooooo.pip27_mean,ooooo.pip28_mean,ooooo.pip29_mean,ooooo.pip30_mean);  
        end
    end
end

for a=1:length(TPHPOS_BNST_L)
    for b=2%:length(names) 
        for c=1:length(pCSname) 
              ooooo=LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.prepip1_mean,ooooo.prepip2_mean,ooooo.prepip3_mean,ooooo.prepip4_mean,ooooo.prepip5_mean,ooooo.prepip6_mean,ooooo.prepip7_mean,ooooo.prepip8_mean,ooooo.prepip9_mean,ooooo.prepip10_mean,ooooo.prepip11_mean,ooooo.prepip12_mean,ooooo.prepip13_mean,ooooo.prepip14_mean,ooooo.prepip15_mean,ooooo.prepip16_mean,ooooo.prepip17_mean,ooooo.prepip18_mean,ooooo.prepip19_mean,ooooo.prepip20_mean,ooooo.prepip21_mean,ooooo.prepip22_mean,ooooo.prepip23_mean,ooooo.prepip24_mean,ooooo.prepip25_mean,ooooo.prepip26_mean,ooooo.prepip27_mean,ooooo.prepip28_mean,ooooo.prepip29_mean,ooooo.prepip30_mean);  
        end
    end
end

%TPHNEG
%CeA

for a=1:length(TPHNEG_CeA_L)
    for b=1%:length(names)
        for c=1:numCS            
            ooooo=LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).altogether=cat(2,ooooo.pip1,ooooo.pip2,ooooo.pip3,ooooo.pip4,ooooo.pip5,ooooo.pip6,ooooo.pip7,ooooo.pip8,ooooo.pip9,ooooo.pip10,ooooo.pip11,ooooo.pip12,ooooo.pip13,ooooo.pip14,ooooo.pip15,ooooo.pip16,ooooo.pip17,ooooo.pip18,ooooo.pip19,ooooo.pip20,ooooo.pip21,ooooo.pip22,ooooo.pip23,ooooo.pip24,ooooo.pip25,ooooo.pip26,ooooo.pip27,ooooo.pip28,ooooo.pip29,ooooo.pip30);
        
            ooooo=LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).altogether=cat(2,ooooo.prepip1,ooooo.prepip2,ooooo.prepip3,ooooo.prepip4,ooooo.prepip5,ooooo.prepip6,ooooo.prepip7,ooooo.prepip8,ooooo.prepip9,ooooo.prepip10,ooooo.prepip11,ooooo.prepip12,ooooo.prepip13,ooooo.prepip14,ooooo.prepip15,ooooo.prepip16,ooooo.prepip17,ooooo.prepip18,ooooo.prepip19,ooooo.prepip20,ooooo.prepip21,ooooo.prepip22,ooooo.prepip23,ooooo.prepip24,ooooo.prepip25,ooooo.prepip26,ooooo.prepip27,ooooo.prepip28,ooooo.prepip29,ooooo.prepip30);
        end
    end
end

for a=1:length(TPHNEG_CeA_L)
    for b=1%:length(names) 
    ooooo=LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime;
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether,ooooo.CS5.altogether,ooooo.CS6.altogether,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs_mean=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs,2,'omitnan');
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4_mean=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4,2,'omitnan');
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10=cat(2,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10_mean=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10,2,'omitnan');
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2_mean=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10.altogether,2);
    end
end


for a=1:length(TPHNEG_CeA_L)
    for b=1%:length(names) 
    ooooo=LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime;
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs=cat(2,ooooo.pCS1.altogether,ooooo.pCS2.altogether,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs_mean=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs,2,'omitnan');
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10=cat(2,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10_mean=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10,2,'omitnan');     
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10only=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10.altogether,2);
    end
end

pip_mean_name={'pip1_mean','pip2_mean','pip3_mean','pip4_mean','pip5_mean','pip6_mean','pip7_mean','pip8_mean','pip9_mean','pip10_mean','pip11_mean','pip12_mean','pip13_mean','pip14_mean','pip15_mean','pip16_mean','pip17_mean','pip18_mean','pip19_mean','pip20_mean','pip21_mean','pip22_mean','pip23_mean','pip24_mean','pip25_mean','pip26_mean','pip27_mean','pip28_mean','pip29_mean','pip30_mean'};
prepip_mean_name={'prepip1_mean','prepip2_mean','prepip3_mean','prepip4_mean','prepip5_mean','prepip6_mean','prepip7_mean','prepip8_mean','prepip9_mean','prepip10_mean','prepip11_mean','prepip12_mean','prepip13_mean','prepip14_mean','prepip15_mean','prepip16_mean','prepip17_mean','prepip18_mean','prepip19_mean','prepip20_mean','prepip21_mean','prepip22_mean','prepip23_mean','prepip24_mean','prepip25_mean','prepip26_mean','prepip27_mean','prepip28_mean','prepip29_mean','prepip30_mean'};

for a=1:length(TPHNEG_CeA_L)
    for b=1%:length(names) 
        for c=1:length(CSname) 
            for d=1:length(Pipname) 
        LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(pip_mean_name{d})=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d}),2,'omitnan');
         end
    end
    end
end


for a=1:length(TPHNEG_CeA_L)
    for b=1%:length(names) 
        for c=1:length(pCSname) 
            for d=1:length(Prepipname) 
        LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(prepip_mean_name{d})=mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d}),2,'omitnan');
         end
    end
    end
end

for a=1:length(TPHNEG_CeA_L)
    for b=1%:length(names) 
        for c=1:length(CSname) 
            ooooo=LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.pip1_mean,ooooo.pip2_mean,ooooo.pip3_mean,ooooo.pip4_mean,ooooo.pip5_mean,ooooo.pip6_mean,ooooo.pip7_mean,ooooo.pip8_mean,ooooo.pip9_mean,ooooo.pip10_mean,ooooo.pip11_mean,ooooo.pip12_mean,ooooo.pip13_mean,ooooo.pip14_mean,ooooo.pip15_mean,ooooo.pip16_mean,ooooo.pip17_mean,ooooo.pip18_mean,ooooo.pip19_mean,ooooo.pip20_mean,ooooo.pip21_mean,ooooo.pip22_mean,ooooo.pip23_mean,ooooo.pip24_mean,ooooo.pip25_mean,ooooo.pip26_mean,ooooo.pip27_mean,ooooo.pip28_mean,ooooo.pip29_mean,ooooo.pip30_mean);  
        end
    end
end

for a=1:length(TPHNEG_CeA_L)
    for b=1%:length(names) 
        for c=1:length(pCSname) 
              ooooo=LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.prepip1_mean,ooooo.prepip2_mean,ooooo.prepip3_mean,ooooo.prepip4_mean,ooooo.prepip5_mean,ooooo.prepip6_mean,ooooo.prepip7_mean,ooooo.prepip8_mean,ooooo.prepip9_mean,ooooo.prepip10_mean,ooooo.prepip11_mean,ooooo.prepip12_mean,ooooo.prepip13_mean,ooooo.prepip14_mean,ooooo.prepip15_mean,ooooo.prepip16_mean,ooooo.prepip17_mean,ooooo.prepip18_mean,ooooo.prepip19_mean,ooooo.prepip20_mean,ooooo.prepip21_mean,ooooo.prepip22_mean,ooooo.prepip23_mean,ooooo.prepip24_mean,ooooo.prepip25_mean,ooooo.prepip26_mean,ooooo.prepip27_mean,ooooo.prepip28_mean,ooooo.prepip29_mean,ooooo.prepip30_mean);  
        end
    end
end

%TPHPOS
%CeA

for a=1:length(TPHPOS_CeA_L)
    for b=1%:length(names)
        for c=1:numCS            
            ooooo=LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).altogether=cat(2,ooooo.pip1,ooooo.pip2,ooooo.pip3,ooooo.pip4,ooooo.pip5,ooooo.pip6,ooooo.pip7,ooooo.pip8,ooooo.pip9,ooooo.pip10,ooooo.pip11,ooooo.pip12,ooooo.pip13,ooooo.pip14,ooooo.pip15,ooooo.pip16,ooooo.pip17,ooooo.pip18,ooooo.pip19,ooooo.pip20,ooooo.pip21,ooooo.pip22,ooooo.pip23,ooooo.pip24,ooooo.pip25,ooooo.pip26,ooooo.pip27,ooooo.pip28,ooooo.pip29,ooooo.pip30);
        
            ooooo=LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).altogether=cat(2,ooooo.prepip1,ooooo.prepip2,ooooo.prepip3,ooooo.prepip4,ooooo.prepip5,ooooo.prepip6,ooooo.prepip7,ooooo.prepip8,ooooo.prepip9,ooooo.prepip10,ooooo.prepip11,ooooo.prepip12,ooooo.prepip13,ooooo.prepip14,ooooo.prepip15,ooooo.prepip16,ooooo.prepip17,ooooo.prepip18,ooooo.prepip19,ooooo.prepip20,ooooo.prepip21,ooooo.prepip22,ooooo.prepip23,ooooo.prepip24,ooooo.prepip25,ooooo.prepip26,ooooo.prepip27,ooooo.prepip28,ooooo.prepip29,ooooo.prepip30);
        end
    end
end

for a=1:length(TPHPOS_CeA_L)
    for b=1%:length(names) 
    ooooo=LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime;
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether,ooooo.CS5.altogether,ooooo.CS6.altogether,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs_mean=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs,2,'omitnan');
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4_mean=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4,2,'omitnan');
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10=cat(2,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10_mean=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10,2,'omitnan');
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2_mean=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10.altogether,2);
    end
end


for a=1:length(TPHPOS_CeA_L)
    for b=1%:length(names) 
    ooooo=LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime;
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs=cat(2,ooooo.pCS1.altogether,ooooo.pCS2.altogether,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs_mean=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs,2,'omitnan');
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10=cat(2,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10_mean=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10,2,'omitnan');     
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10only=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10.altogether,2);
    end
end

pip_mean_name={'pip1_mean','pip2_mean','pip3_mean','pip4_mean','pip5_mean','pip6_mean','pip7_mean','pip8_mean','pip9_mean','pip10_mean','pip11_mean','pip12_mean','pip13_mean','pip14_mean','pip15_mean','pip16_mean','pip17_mean','pip18_mean','pip19_mean','pip20_mean','pip21_mean','pip22_mean','pip23_mean','pip24_mean','pip25_mean','pip26_mean','pip27_mean','pip28_mean','pip29_mean','pip30_mean'};
prepip_mean_name={'prepip1_mean','prepip2_mean','prepip3_mean','prepip4_mean','prepip5_mean','prepip6_mean','prepip7_mean','prepip8_mean','prepip9_mean','prepip10_mean','prepip11_mean','prepip12_mean','prepip13_mean','prepip14_mean','prepip15_mean','prepip16_mean','prepip17_mean','prepip18_mean','prepip19_mean','prepip20_mean','prepip21_mean','prepip22_mean','prepip23_mean','prepip24_mean','prepip25_mean','prepip26_mean','prepip27_mean','prepip28_mean','prepip29_mean','prepip30_mean'};

for a=1:length(TPHPOS_CeA_L)
    for b=1%:length(names) 
        for c=1:length(CSname) 
            for d=1:length(Pipname) 
        LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(pip_mean_name{d})=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d}),2,'omitnan');
         end
    end
    end
end


for a=1:length(TPHPOS_CeA_L)
    for b=1%:length(names) 
        for c=1:length(pCSname) 
            for d=1:length(Prepipname) 
        LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(prepip_mean_name{d})=mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d}),2,'omitnan');
         end
    end
    end
end

for a=1:length(TPHPOS_CeA_L)
    for b=1%:length(names) 
        for c=1:length(CSname) 
            ooooo=LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.pip1_mean,ooooo.pip2_mean,ooooo.pip3_mean,ooooo.pip4_mean,ooooo.pip5_mean,ooooo.pip6_mean,ooooo.pip7_mean,ooooo.pip8_mean,ooooo.pip9_mean,ooooo.pip10_mean,ooooo.pip11_mean,ooooo.pip12_mean,ooooo.pip13_mean,ooooo.pip14_mean,ooooo.pip15_mean,ooooo.pip16_mean,ooooo.pip17_mean,ooooo.pip18_mean,ooooo.pip19_mean,ooooo.pip20_mean,ooooo.pip21_mean,ooooo.pip22_mean,ooooo.pip23_mean,ooooo.pip24_mean,ooooo.pip25_mean,ooooo.pip26_mean,ooooo.pip27_mean,ooooo.pip28_mean,ooooo.pip29_mean,ooooo.pip30_mean);  
        end
    end
end

for a=1:length(TPHPOS_CeA_L)
    for b=1%:length(names) 
        for c=1:length(pCSname) 
              ooooo=LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.prepip1_mean,ooooo.prepip2_mean,ooooo.prepip3_mean,ooooo.prepip4_mean,ooooo.prepip5_mean,ooooo.prepip6_mean,ooooo.prepip7_mean,ooooo.prepip8_mean,ooooo.prepip9_mean,ooooo.prepip10_mean,ooooo.prepip11_mean,ooooo.prepip12_mean,ooooo.prepip13_mean,ooooo.prepip14_mean,ooooo.prepip15_mean,ooooo.prepip16_mean,ooooo.prepip17_mean,ooooo.prepip18_mean,ooooo.prepip19_mean,ooooo.prepip20_mean,ooooo.prepip21_mean,ooooo.prepip22_mean,ooooo.prepip23_mean,ooooo.prepip24_mean,ooooo.prepip25_mean,ooooo.prepip26_mean,ooooo.prepip27_mean,ooooo.prepip28_mean,ooooo.prepip29_mean,ooooo.prepip30_mean);  
        end
    end
end

%% plot the pwr graph for each animal
    %bnst
for a=1:length(TPHNEG_BNST_L)
figure (a)
    for b=2 
     subplot(1,2,(b-1))
        for c=1:10
           
uuuu=LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
plot(mean(uuuu.prepip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHNEG ' char(TPHNEG_BNST_L(a)), 'all pCSs pips ' char(names(b))])


   subplot(1,2,b)
        for c=1:10
           
uuuu=LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
plot(mean(uuuu.pip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHNEG ' char(TPHNEG_BNST_L(a)), 'all CSs pips ' char(names(b))])
        end
end

for a=1:length(TPHPOS_BNST_L)
figure (a+40)
    for b=2 
     subplot(1,2,(b-1))
        for c=1:10
           
uuuu=LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
plot(mean(uuuu.prepip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHPOS ' char(TPHNEG_BNST_L(a)), 'all pCSs pips ' char(names(b))])


   subplot(1,2,b)
        for c=1:10
           
uuuu=LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
plot(mean(uuuu.pip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHPOS ' char(TPHPOS_BNST_L(a)), 'all CSs pips ' char(names(b))])
        end
    end

%% plot the pwr graph for each animal
    %cea
for a=1:length(TPHNEG_CeA_L)
figure (a)
    for b=1 
     subplot(1,2,(b))
        for c=1:10
           
uuuu=LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
plot(mean(uuuu.prepip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHNEG ' char(TPHNEG_CeA_L(a)), 'all pCSs pips ' char(names(b))])


   subplot(1,2,b+1)
        for c=1:10
           
uuuu=LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
plot(mean(uuuu.pip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHNEG ' char(TPHNEG_CeA_L(a)), 'all CSs pips ' char(names(b))])
        end
end

for a=1:length(TPHPOS_CeA_L)
figure (a+10)
    for b=1 
     subplot(1,2,(b))
        for c=1:10
           
uuuu=LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
plot(mean(uuuu.prepip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHPOS ' char(TPHPOS_CeA_L(a)), 'all pCSs pips ' char(names(b))])


   subplot(1,2,b+1)
        for c=1:10
           
uuuu=LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
plot(mean(uuuu.pip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHPOS ' char(TPHPOS_CeA_L(a)), 'all CSs pips ' char(names(b))])
        end
    end

 %% pwr calcs
%cea
for a=1:length(TPHNEG_CeA_L)
     for n=1%:length(names)
    fooTPHNEG_CeA=LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{n}).pips.Pips_freqtime;
Pwr.TPHNEG.(names{n}).CS1_4_mean(:,a)=fooTPHNEG_CeA.CS1_4_mean;
Pwr.TPHNEG.(names{n}).CS7_10_mean(:,a)=fooTPHNEG_CeA.CS7_10_mean;
Pwr.TPHNEG.(names{n}).all_CSs_mean(:,a)=fooTPHNEG_CeA.all_CSs_mean;
Pwr.TPHNEG.(names{n}).CS1_mean(:,a)=fooTPHNEG_CeA.CS1only;
Pwr.TPHNEG.(names{n}).CS2_mean(:,a)=fooTPHNEG_CeA.CS2only;
Pwr.TPHNEG.(names{n}).CS3_mean(:,a)=fooTPHNEG_CeA.CS3only;
Pwr.TPHNEG.(names{n}).CS4_mean(:,a)=fooTPHNEG_CeA.CS4only;
Pwr.TPHNEG.(names{n}).CS5_mean(:,a)=fooTPHNEG_CeA.CS5only;
Pwr.TPHNEG.(names{n}).CS6_mean(:,a)=fooTPHNEG_CeA.CS6only;
Pwr.TPHNEG.(names{n}).CS7_mean(:,a)=fooTPHNEG_CeA.CS7only;
Pwr.TPHNEG.(names{n}).CS8_mean(:,a)=fooTPHNEG_CeA.CS8only;
Pwr.TPHNEG.(names{n}).CS9_mean(:,a)=fooTPHNEG_CeA.CS9only;
Pwr.TPHNEG.(names{n}).CS10_mean(:,a)=fooTPHNEG_CeA.CS10only;

    fooTPHNEG_CeA=LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{n}).pips.Prepips_freqtime; 
Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,a)=fooTPHNEG_CeA.all_pCSs_mean;

Pwr.TPHNEG.(names{n}).normCS1_4_mean=Pwr.TPHNEG.(names{n}).CS1_4_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).normCS7_10_mean=Pwr.TPHNEG.(names{n}).CS7_10_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:); 
Pwr.TPHNEG.(names{n}).norm_all_CSs_mean=Pwr.TPHNEG.(names{n}).all_CSs_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS1_mean=Pwr.TPHNEG.(names{n}).CS1_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS2_mean=Pwr.TPHNEG.(names{n}).CS2_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS3_mean=Pwr.TPHNEG.(names{n}).CS3_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS4_mean=Pwr.TPHNEG.(names{n}).CS4_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS5_mean=Pwr.TPHNEG.(names{n}).CS5_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS6_mean=Pwr.TPHNEG.(names{n}).CS6_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS7_mean=Pwr.TPHNEG.(names{n}).CS7_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS8_mean=Pwr.TPHNEG.(names{n}).CS8_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS9_mean=Pwr.TPHNEG.(names{n}).CS9_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS10_mean=Pwr.TPHNEG.(names{n}).CS10_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
     end
end

for a=1:length(TPHPOS_CeA_L)
     for n=1%:length(names)
    fooTPHPOS_CeA=LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{n}).pips.Pips_freqtime;
Pwr.TPHPOS.(names{n}).CS1_4_mean(:,a)=fooTPHPOS_CeA.CS1_4_mean;
Pwr.TPHPOS.(names{n}).CS7_10_mean(:,a)=fooTPHPOS_CeA.CS7_10_mean;
Pwr.TPHPOS.(names{n}).all_CSs_mean(:,a)=fooTPHPOS_CeA.all_CSs_mean;
Pwr.TPHPOS.(names{n}).CS1_mean(:,a)=fooTPHPOS_CeA.CS1only;
Pwr.TPHPOS.(names{n}).CS2_mean(:,a)=fooTPHPOS_CeA.CS2only;
Pwr.TPHPOS.(names{n}).CS3_mean(:,a)=fooTPHPOS_CeA.CS3only;
Pwr.TPHPOS.(names{n}).CS4_mean(:,a)=fooTPHPOS_CeA.CS4only;
Pwr.TPHPOS.(names{n}).CS5_mean(:,a)=fooTPHPOS_CeA.CS5only;
Pwr.TPHPOS.(names{n}).CS6_mean(:,a)=fooTPHPOS_CeA.CS6only;
Pwr.TPHPOS.(names{n}).CS7_mean(:,a)=fooTPHPOS_CeA.CS7only;
Pwr.TPHPOS.(names{n}).CS8_mean(:,a)=fooTPHPOS_CeA.CS8only;
Pwr.TPHPOS.(names{n}).CS9_mean(:,a)=fooTPHPOS_CeA.CS9only;
Pwr.TPHPOS.(names{n}).CS10_mean(:,a)=fooTPHPOS_CeA.CS10only;

    fooTPHPOS_CeA=LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{n}).pips.Prepips_freqtime; 
Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,a)=fooTPHPOS_CeA.all_pCSs_mean;

Pwr.TPHPOS.(names{n}).normCS1_4_mean=Pwr.TPHPOS.(names{n}).CS1_4_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).normCS7_10_mean=Pwr.TPHPOS.(names{n}).CS7_10_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:); 
Pwr.TPHPOS.(names{n}).norm_all_CSs_mean=Pwr.TPHPOS.(names{n}).all_CSs_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS1_mean=Pwr.TPHPOS.(names{n}).CS1_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS2_mean=Pwr.TPHPOS.(names{n}).CS2_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS3_mean=Pwr.TPHPOS.(names{n}).CS3_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS4_mean=Pwr.TPHPOS.(names{n}).CS4_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS5_mean=Pwr.TPHPOS.(names{n}).CS5_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS6_mean=Pwr.TPHPOS.(names{n}).CS6_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS7_mean=Pwr.TPHPOS.(names{n}).CS7_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS8_mean=Pwr.TPHPOS.(names{n}).CS8_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS9_mean=Pwr.TPHPOS.(names{n}).CS9_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS10_mean=Pwr.TPHPOS.(names{n}).CS10_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);

end
end

 %bnst
fooTPHNEG_BNST=LFPsall.TPHNEG;
for a=1:length(TPHNEG_BNST_L)
     for n=2%:length(names)
    fooTPHNEG_BNST=LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{n}).pips.Pips_freqtime;
Pwr.TPHNEG.(names{n}).CS1_4_mean(:,a)=fooTPHNEG_BNST.CS1_4_mean;
Pwr.TPHNEG.(names{n}).CS7_10_mean(:,a)=fooTPHNEG_BNST.CS7_10_mean;
Pwr.TPHNEG.(names{n}).all_CSs_mean(:,a)=fooTPHNEG_BNST.all_CSs_mean;
Pwr.TPHNEG.(names{n}).CS1_mean(:,a)=fooTPHNEG_BNST.CS1only;
Pwr.TPHNEG.(names{n}).CS2_mean(:,a)=fooTPHNEG_BNST.CS2only;
Pwr.TPHNEG.(names{n}).CS3_mean(:,a)=fooTPHNEG_BNST.CS3only;
Pwr.TPHNEG.(names{n}).CS4_mean(:,a)=fooTPHNEG_BNST.CS4only;
Pwr.TPHNEG.(names{n}).CS5_mean(:,a)=fooTPHNEG_BNST.CS5only;
Pwr.TPHNEG.(names{n}).CS6_mean(:,a)=fooTPHNEG_BNST.CS6only;
Pwr.TPHNEG.(names{n}).CS7_mean(:,a)=fooTPHNEG_BNST.CS7only;
Pwr.TPHNEG.(names{n}).CS8_mean(:,a)=fooTPHNEG_BNST.CS8only;
Pwr.TPHNEG.(names{n}).CS9_mean(:,a)=fooTPHNEG_BNST.CS9only;
Pwr.TPHNEG.(names{n}).CS10_mean(:,a)=fooTPHNEG_BNST.CS10only;

    fooTPHNEG_BNST=LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{n}).pips.Prepips_freqtime; 
Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,a)=fooTPHNEG_BNST.all_pCSs_mean;

Pwr.TPHNEG.(names{n}).normCS1_4_mean=Pwr.TPHNEG.(names{n}).CS1_4_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).normCS7_10_mean=Pwr.TPHNEG.(names{n}).CS7_10_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:); 
Pwr.TPHNEG.(names{n}).norm_all_CSs_mean=Pwr.TPHNEG.(names{n}).all_CSs_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS1_mean=Pwr.TPHNEG.(names{n}).CS1_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS2_mean=Pwr.TPHNEG.(names{n}).CS2_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS3_mean=Pwr.TPHNEG.(names{n}).CS3_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS4_mean=Pwr.TPHNEG.(names{n}).CS4_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS5_mean=Pwr.TPHNEG.(names{n}).CS5_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS6_mean=Pwr.TPHNEG.(names{n}).CS6_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS7_mean=Pwr.TPHNEG.(names{n}).CS7_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS8_mean=Pwr.TPHNEG.(names{n}).CS8_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS9_mean=Pwr.TPHNEG.(names{n}).CS9_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHNEG.(names{n}).norm_CS10_mean=Pwr.TPHNEG.(names{n}).CS10_mean(:,:)./Pwr.TPHNEG.(names{n}).all_pCSs_mean(:,:);

end
end

fooTPHPOS_BNST=LFPsall.TPHPOS;
for a=1:length(TPHPOS_BNST_L)
     for n=2%:length(names)
    fooTPHPOS_BNST=LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{n}).pips.Pips_freqtime;
Pwr.TPHPOS.(names{n}).CS1_4_mean(:,a)=fooTPHPOS_BNST.CS1_4_mean;
Pwr.TPHPOS.(names{n}).CS7_10_mean(:,a)=fooTPHPOS_BNST.CS7_10_mean;
Pwr.TPHPOS.(names{n}).all_CSs_mean(:,a)=fooTPHPOS_BNST.all_CSs_mean;
Pwr.TPHPOS.(names{n}).CS1_mean(:,a)=fooTPHPOS_BNST.CS1only;
Pwr.TPHPOS.(names{n}).CS2_mean(:,a)=fooTPHPOS_BNST.CS2only;
Pwr.TPHPOS.(names{n}).CS3_mean(:,a)=fooTPHPOS_BNST.CS3only;
Pwr.TPHPOS.(names{n}).CS4_mean(:,a)=fooTPHPOS_BNST.CS4only;
Pwr.TPHPOS.(names{n}).CS5_mean(:,a)=fooTPHPOS_BNST.CS5only;
Pwr.TPHPOS.(names{n}).CS6_mean(:,a)=fooTPHPOS_BNST.CS6only;
Pwr.TPHPOS.(names{n}).CS7_mean(:,a)=fooTPHPOS_BNST.CS7only;
Pwr.TPHPOS.(names{n}).CS8_mean(:,a)=fooTPHPOS_BNST.CS8only;
Pwr.TPHPOS.(names{n}).CS9_mean(:,a)=fooTPHPOS_BNST.CS9only;
Pwr.TPHPOS.(names{n}).CS10_mean(:,a)=fooTPHPOS_BNST.CS10only;

    fooTPHPOS_BNST=LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{n}).pips.Prepips_freqtime; 
Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,a)=fooTPHPOS_BNST.all_pCSs_mean;

Pwr.TPHPOS.(names{n}).normCS1_4_mean=Pwr.TPHPOS.(names{n}).CS1_4_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).normCS7_10_mean=Pwr.TPHPOS.(names{n}).CS7_10_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:); 
Pwr.TPHPOS.(names{n}).norm_all_CSs_mean=Pwr.TPHPOS.(names{n}).all_CSs_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS1_mean=Pwr.TPHPOS.(names{n}).CS1_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS2_mean=Pwr.TPHPOS.(names{n}).CS2_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS3_mean=Pwr.TPHPOS.(names{n}).CS3_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS4_mean=Pwr.TPHPOS.(names{n}).CS4_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS5_mean=Pwr.TPHPOS.(names{n}).CS5_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS6_mean=Pwr.TPHPOS.(names{n}).CS6_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS7_mean=Pwr.TPHPOS.(names{n}).CS7_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS8_mean=Pwr.TPHPOS.(names{n}).CS8_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS9_mean=Pwr.TPHPOS.(names{n}).CS9_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwr.TPHPOS.(names{n}).norm_CS10_mean=Pwr.TPHPOS.(names{n}).CS10_mean(:,:)./Pwr.TPHPOS.(names{n}).all_pCSs_mean(:,:);

end
end

%%
what={'all_CSs_mean', 'all_pCSs_mean','CS1_4_mean','CS7_10_mean'};

for n=1:length(names)
    for t=1:length(what)
     Pwr_all.(names{n}).TPHPOS.(what{t})=(Pwr.TPHPOS.(names{n}).(what{t}));
     Pwr_all.(names{n}).TPHNEG.(what{t})=(Pwr.TPHNEG.(names{n}).(what{t}));
    end
end

Freq=Freq_pw;
gamma_freq_idx_90_140=find(Freq>=90 & Freq<=140);

figure
n=1 % cea

%tph- ACSF
 subplot(3,4,1)     
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{1})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]); ylim([1 80000]) ;title([char(names(n)), ' RECALL TPH-' ]); 
        
        subplot(3,4,5)
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{3})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 80000]); title([char(names(n)), ' RECALL TPH-' ]);
        
        subplot(3,4,9)
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{4})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 80000]); title([char(names(n)), ' RECALL TPH-']); 

%TPH+ ACSF
     subplot(3,4,2)     
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{1})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]); ylim([1 80000]) ;title([char(names(n)), ' RECALL TPH+' ]); 
        
        subplot(3,4,6)
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{3})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 80000]); title([char(names(n)), ' RECALL TPH+' ]);
        
        subplot(3,4,10)
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{4})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 80000]); title([char(names(n)), ' RECALL TPH+']); 
       

        %TPH- ACSF
        fooTPHNEG=struct2cell(Pwr_all.(names{n}).TPHNEG);
        %CeA L
        subplot(3,4,3)
        CS_SEM = std(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH- all CSs ', ' ','TPH- all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,7)
        CS_SEM = std(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
         ylim([0 80000]);legend('TPH- CS1-4 ', ' ','TPH- all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,11)
        CS_SEM = std(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH- CS7-10 ', ' ','TPH- all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
   
        %TPH+ acsf
       fooTPHPOS=struct2cell(Pwr_all.(names{n}).TPHPOS);
        %CeA L
        subplot(3,4,4)
        CS_SEM = std(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH+ all CSs ', ' ','TPH+ all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,8)
        CS_SEM = std(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
         ylim([0 80000]);legend('TPH+ CS1-4 ', ' ','TPH+ all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,12)
        CS_SEM = std(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH+ CS7-10 ', ' ','TPH+ all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
     

%% 
%Absolute Power for each region, CS and pCS
what={'all_CSs_mean', 'all_pCSs_mean','CS1_4_mean','CS7_10_mean'};

Freq=Freq_pw;
gamma_freq_idx_90_140=find(Freq>=70 & Freq<=150);

figure
n=2 % BNST

%tph- ACSF
 subplot(3,4,1)     
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{1})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]); ylim([1 80000]) ;title([char(names(n)), ' RECALL TPH-' ]); 
        
        subplot(3,4,5)
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{3})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 80000]); title([char(names(n)), ' RECALL TPH-' ]);
        
        subplot(3,4,9)
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{4})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 80000]); title([char(names(n)), ' RECALL TPH-']); 

%TPH+ ACSF
     subplot(3,4,2)     
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{1})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]); ylim([1 80000]) ;title([char(names(n)), ' RECALL TPH+' ]); 
        
        subplot(3,4,6)
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{3})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 80000]); title([char(names(n)), ' RECALL TPH+' ]);
        
        subplot(3,4,10)
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{4})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 80000]); title([char(names(n)), ' RECALL TPH+']); 
       

        %TPH- ACSF
        fooTPHNEG=struct2cell(Pwr_all.(names{n}).TPHNEG);
        %CeA L
        subplot(3,4,3)
        CS_SEM = std(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH- all CSs ', ' ','TPH- all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,7)
        CS_SEM = std(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
         ylim([0 80000]);legend('TPH- CS1-4 ', ' ','TPH- all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,11)
        CS_SEM = std(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH- CS7-10 ', ' ','TPH- all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
   
        %TPH+ acsf
       fooTPHPOS=struct2cell(Pwr_all.(names{n}).TPHPOS);
        %CeA L
        subplot(3,4,4)
        CS_SEM = std(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH+ all CSs ', ' ','TPH+ all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,8)
        CS_SEM = std(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
         ylim([0 80000]);legend('TPH+ CS1-4 ', ' ','TPH+ all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,12)
        CS_SEM = std(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH+ CS7-10 ', ' ','TPH+ all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
           
%%
%NORM Power for each region, CS and pCS
clear what
clear Pwr_all
what={'norm_all_CSs_mean', 'normCS1_4_mean', 'normCS7_10_mean'};

for n=1:length(names)
    for t=1:length(what)
     Pwr_all.(names{n}).TPHPOS.(what{t})=(Pwr.TPHPOS.(names{n}).(what{t}));
     Pwr_all.(names{n}).TPHNEG.(what{t})=(Pwr.TPHNEG.(names{n}).(what{t}));
    end
end

figure(23)   
n=1 ;%CeA_L
        subplot(3,4,1)     
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{1})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{1})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')
       
        subplot(3,4,5)  
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{2})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{2})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')
        
        subplot(3,4,9)
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{3})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{3})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')

subplot(3,4,2)
        fooTPHNEG=struct2cell(Pwr_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')  

        fooTPHPOS=struct2cell(Pwr_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1.5]);legend('TPH- all CS ',' ','TPH+ all CS ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
      
        subplot(3,4,6)
        fooTPHNEG=struct2cell(Pwr_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
  
        fooTPHPOS=struct2cell(Pwr_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
           ylim([0 1.5]);legend('TPH- CS1-4 ',' ','TPH+ CS1-4 ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
      
        subplot(3,4,10)
        fooTPHNEG=struct2cell(Pwr_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        

        fooTPHPOS=struct2cell(Pwr_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1.5]);legend('TPH- CS7-10 ', ' ','TPH+ CS7-10 ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])       

        n=2 ;%bnst

        subplot(3,4,3)     
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{1})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{1})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')
       
        subplot(3,4,7)  
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{2})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{2})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')
        
        subplot(3,4,11)
        errorbarplot_joe(Freq_pw,Pwr.TPHNEG.(names{n}).(what{3})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwr.TPHPOS.(names{n}).(what{3})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')

subplot(3,4,4)
        fooTPHNEG=struct2cell(Pwr_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')  

        fooTPHPOS=struct2cell(Pwr_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1.5]);legend('TPH- all CS ',' ','TPH+ all CS ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
      
        subplot(3,4,8)
        fooTPHNEG=struct2cell(Pwr_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
  
        fooTPHPOS=struct2cell(Pwr_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
           ylim([0 1.5]);legend('TPH- CS1-4 ',' ','TPH+ CS1-4 ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
      
        subplot(3,4,12)
        fooTPHNEG=struct2cell(Pwr_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        
        fooTPHPOS=struct2cell(Pwr_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1.5]);legend('TPH- CS7-10 ', ' ','TPH+ CS7-10 ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])       

%% NaN
%MANUALLY DUPLICATE LFPSALL AND RENAME TO LFPSALLNAN

%cea
forcutoff1=struct('TPH645_RECALL',[],'TPH646_RECALL',[],'TPH647_RECALL',[],'TPH891_RECALL',[],'TPH1006_RECALL',[],'TPH1406_RECALL',[],'TPH1509_RECALL',[]);

for a=1:length(TPHNEG_CeA_L)
      forcutoff1.(TPHNEG_CeA_L{a})=struct('CeA',[]);

end

for a=1:length(TPHNEG_CeA_L)
        for b=1%1:length(names)
         forcutoff1.(TPHNEG_CeA_L{a}).(names{b})=cat(2,LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10.meaneachpip_altogether(73:154,:))
    end
end

for a=1:length(TPHNEG_CeA_L)
        for b=1%1:length(names)
        cutoff2.(TPHNEG_CeA_L{a}).(names{b})=median(forcutoff1.(TPHNEG_CeA_L{a}).(names{b}),2)
    end
end

for a=1:length(TPHNEG_CeA_L)
        for b=1%1:length(names)
        cutoff3.(TPHNEG_CeA_L{a}).(names{b})=(mean(cutoff2.(TPHNEG_CeA_L{a}).(names{b})*15));
    end
end

 for a=1:length(TPHNEG_CeA_L)
     for b=1%1:length(names)
          for c=1:length(CSname)
              for d=1:length(Pipname)

 	xx_ind1=mean(mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d})(73:154,:),2))>(cutoff3.(TPHNEG_CeA_L{a}).(names{b}));
    find(xx_ind1==1)
    if find(xx_ind1==1) ~= 0
     
   LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d})(:,:)=nan;       
    end
              end
          end
     end
 end

 %preCS

forcutoff1_pre=struct('TPH645_RECALL',[],'TPH646_RECALL',[],'TPH647_RECALL',[],'TPH891_RECALL',[],'TPH1006_RECALL',[],'TPH1406_RECALL',[],'TPH1509_RECALL',[]);

for a=1:length(TPHNEG_CeA_L)
      forcutoff1_pre.(TPHNEG_CeA_L{a})=struct('CeA',[]);

end

for a=1:length(TPHNEG_CeA_L)
        for b=1%1:length(names)
         forcutoff1_pre.(TPHNEG_CeA_L{a}).(names{b})=cat(2,LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10.meaneachpip_altogether(73:154,:))
    end
end

for a=1:length(TPHNEG_CeA_L)
        for b=1%1:length(names)
        cutoff2_pre.(TPHNEG_CeA_L{a}).(names{b})=median(forcutoff1_pre.(TPHNEG_CeA_L{a}).(names{b}),2)
    end
end

for a=1:length(TPHNEG_CeA_L)
        for b=1%1:length(names)
        cutoff3_pre.(TPHNEG_CeA_L{a}).(names{b})=(mean(cutoff2_pre.(TPHNEG_CeA_L{a}).(names{b})*15));
    end
end

 for a=1:length(TPHNEG_CeA_L)
     for b=1%1:length(names)
          for c=1:length(pCSname)
              for d=1:length(Prepipname)

 	xx_ind1=mean(mean(LFPsall.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d})(73:154,:),2))>(cutoff3_pre.(TPHNEG_CeA_L{a}).(names{b}));
    find(xx_ind1==1)
    if find(xx_ind1==1) ~= 0
     
   LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d})(:,:)=nan;       
    end
              end
          end
     end
 end

clear forcuttoff1
clear forcuttoff2
clear forcuttoff3
clear forcuttoff1_pre
clear forcuttoff2_pre
clear forcuttoff3_pre

 %TPH+
forcutoff1=struct('TPH655_RECALL',[],'TPH691_RECALL',[],'TPH905_RECALL',[],'TPH907_RECALL',[],'TPH1008_RECALL',[],'TPH1026_RECALL',[],'TPH1390_RECALL',[]);

for a=1:length(TPHPOS_CeA_L)
      forcutoff1.(TPHPOS_CeA_L{a})=struct('CeA',[]);

end

for a=1:length(TPHPOS_CeA_L)
        for b=1%1:length(names)
         forcutoff1.(TPHPOS_CeA_L{a}).(names{b})=cat(2,LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10.meaneachpip_altogether(73:154,:))
    end
end

for a=1:length(TPHPOS_CeA_L)
        for b=1%1:length(names)
        cutoff2.(TPHPOS_CeA_L{a}).(names{b})=median(forcutoff1.(TPHPOS_CeA_L{a}).(names{b}),2)
    end
end

for a=1:length(TPHPOS_CeA_L)
        for b=1%1:length(names)
        cutoff3.(TPHPOS_CeA_L{a}).(names{b})=(mean(cutoff2.(TPHPOS_CeA_L{a}).(names{b})*15));
    end
end

 for a=1:length(TPHPOS_CeA_L)
     for b=1%1:length(names)
          for c=1:length(CSname)
              for d=1:length(Pipname)

 	xx_ind1=mean(mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d})(73:154,:),2))>(cutoff3.(TPHPOS_CeA_L{a}).(names{b}));
    find(xx_ind1==1)
    if find(xx_ind1==1) ~= 0
     
   LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d})(:,:)=nan;       
    end
              end
          end
     end
 end

 %preCS
forcutoff1_pre=struct('TPH655_RECALL',[],'TPH691_RECALL',[],'TPH905_RECALL',[],'TPH907_RECALL',[],'TPH1008_RECALL',[],'TPH1026_RECALL',[],'TPH1390_RECALL',[]);

for a=1:length(TPHPOS_CeA_L)
      forcutoff1_pre.(TPHPOS_CeA_L{a})=struct('CeA',[]);

end

for a=1:length(TPHPOS_CeA_L)
        for b=1%1:length(names)
         forcutoff1_pre.(TPHPOS_CeA_L{a}).(names{b})=cat(2,LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10.meaneachpip_altogether(73:154,:))
    end
end

for a=1:length(TPHPOS_CeA_L)
        for b=1%1:length(names)
        cutoff2_pre.(TPHPOS_CeA_L{a}).(names{b})=median(forcutoff1_pre.(TPHPOS_CeA_L{a}).(names{b}),2)
    end
end

for a=1:length(TPHPOS_CeA_L)
        for b=1%1:length(names)
        cutoff3_pre.(TPHPOS_CeA_L{a}).(names{b})=(mean(cutoff2_pre.(TPHPOS_CeA_L{a}).(names{b})*15));
    end
end

 for a=1:length(TPHPOS_CeA_L)
     for b=1%1:length(names)
          for c=1:length(pCSname)
              for d=1:length(Prepipname)

 	xx_ind1=mean(mean(LFPsall.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d})(73:154,:),2))>(cutoff3_pre.(TPHPOS_CeA_L{a}).(names{b}));
    find(xx_ind1==1)
    if find(xx_ind1==1) ~= 0
     
   LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d})(:,:)=nan;       
    end
              end
          end
     end
 end

clear forcuttoff1
clear forcuttoff2
clear forcuttoff3
clear forcuttoff1_pre
clear forcuttoff2_pre
clear forcuttoff3_pre

%bnst
forcutoff1=struct('TPH645_RECALL',[],'TPH646_RECALL',[],'TPH647_RECALL',[],'TPH662_RECALL',[],'TPH891_RECALL',[],'TPH906_RECALL',[],'TPH1006_RECALL',[],'TPH1025_RECALL',[], 'TPH1398_RECALL',[],'TPH1406_RECALL',[]);

for a=1:length(TPHNEG_BNST_L)
      forcutoff1.(TPHNEG_BNST_L{a})=struct('BNST',[]);

end

for a=1:length(TPHNEG_BNST_L)
        for b=2%1:length(names)
         forcutoff1.(TPHNEG_BNST_L{a}).(names{b})=cat(2,LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10.meaneachpip_altogether(73:154,:))
    end
end

for a=1:length(TPHNEG_BNST_L)
        for b=2%1:length(names)
        cutoff2.(TPHNEG_BNST_L{a}).(names{b})=median(forcutoff1.(TPHNEG_BNST_L{a}).(names{b}),2)
    end
end

for a=1:length(TPHNEG_BNST_L)
        for b=2%1:length(names)
        cutoff3.(TPHNEG_BNST_L{a}).(names{b})=(mean(cutoff2.(TPHNEG_BNST_L{a}).(names{b})*15));
    end
end

 for a=1:length(TPHNEG_BNST_L)
     for b=2%1:length(names)
          for c=1:length(CSname)
              for d=1:length(Pipname)

 	xx_ind1=mean(mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d})(73:154,:),2))>(cutoff3.(TPHNEG_BNST_L{a}).(names{b}));
    find(xx_ind1==1)
    if find(xx_ind1==1) ~= 0
     
   LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d})(:,:)=nan;       
    end
              end
          end
     end
 end

 %preCS

forcutoff1_pre=struct('TPH645_RECALL',[],'TPH646_RECALL',[],'TPH647_RECALL',[],'TPH662_RECALL',[],'TPH891_RECALL',[],'TPH906_RECALL',[],'TPH1006_RECALL',[],'TPH1025_RECALL',[], 'TPH1398_RECALL',[],'TPH1406_RECALL',[]);

for a=1:length(TPHNEG_BNST_L)
      forcutoff1_pre.(TPHNEG_BNST_L{a})=struct('BNST',[]);

end

for a=1:length(TPHNEG_BNST_L)
        for b=2%1:length(names)
         forcutoff1_pre.(TPHNEG_BNST_L{a}).(names{b})=cat(2,LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9.meaneachpip_altogether(73:154,:),LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10.meaneachpip_altogether(73:154,:))
    end
end

for a=1:length(TPHNEG_BNST_L)
        for b=2%1:length(names)
        cutoff2_pre.(TPHNEG_BNST_L{a}).(names{b})=median(forcutoff1_pre.(TPHNEG_BNST_L{a}).(names{b}),2)
    end
end

for a=1:length(TPHNEG_BNST_L)
        for b=2%1:length(names)
        cutoff3_pre.(TPHNEG_BNST_L{a}).(names{b})=(mean(cutoff2_pre.(TPHNEG_BNST_L{a}).(names{b})*15));
    end
end

 for a=1:length(TPHNEG_BNST_L)
     for b=2%1:length(names)
          for c=1:length(pCSname)
              for d=1:length(Prepipname)

 	xx_ind1=mean(mean(LFPsall.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d})(73:154,:),2))>(cutoff3_pre.(TPHNEG_BNST_L{a}).(names{b}));
    find(xx_ind1==1)
    if find(xx_ind1==1) ~= 0
     
   LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d})(:,:)=nan;       
    end
              end
          end
     end
 end

clear forcuttoff1
clear forcuttoff2
clear forcuttoff3
clear forcuttoff1_pre
clear forcuttoff2_pre
clear forcuttoff3_pre

 %TPH+
forcutoff1=struct('TPH394_RECALL',[],'TPH655_RECALL',[],'TPH657_RECALL',[],'TPH905_RECALL',[],'TPH907_RECALL',[],'TPH1008_RECALL',[],'TPH1026_RECALL',[],'TPH1389_RECALL',[],'TPH1390_RECALL',[],'TPH1397_RECALL',[],'TPH1399_RECALL',[]);

for a=1:length(TPHPOS_BNST_L)
      forcutoff1.(TPHPOS_BNST_L{a})=struct('BNST',[]);

end

for a=1:length(TPHPOS_BNST_L)
        for b=2%1:length(names)
         forcutoff1.(TPHPOS_BNST_L{a}).(names{b})=cat(2,LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10.meaneachpip_altogether(73:154,:))
    end
end

for a=1:length(TPHPOS_BNST_L)
        for b=2%1:length(names)
        cutoff2.(TPHPOS_BNST_L{a}).(names{b})=median(forcutoff1.(TPHPOS_BNST_L{a}).(names{b}),2)
    end
end

for a=1:length(TPHPOS_BNST_L)
        for b=2%1:length(names)
        cutoff3.(TPHPOS_BNST_L{a}).(names{b})=(mean(cutoff2.(TPHPOS_BNST_L{a}).(names{b})*15));
    end
end

 for a=1:length(TPHPOS_BNST_L)
     for b=2%1:length(names)
          for c=1:length(CSname)
              for d=1:length(Pipname)

 	xx_ind1=mean(mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d})(73:154,:),2))>(cutoff3.(TPHPOS_BNST_L{a}).(names{b}));
    find(xx_ind1==1)
    if find(xx_ind1==1) ~= 0
     
   LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d})(:,:)=nan;       
    end
              end
          end
     end
 end

 %preCS
forcutoff1_pre=struct('TPH394_RECALL',[],'TPH655_RECALL',[],'TPH657_RECALL',[],'TPH905_RECALL',[],'TPH907_RECALL',[],'TPH1008_RECALL',[],'TPH1026_RECALL',[],'TPH1389_RECALL',[],'TPH1390_RECALL',[],'TPH1397_RECALL',[],'TPH1399_RECALL',[]);

for a=1:length(TPHPOS_BNST_L)
      forcutoff1_pre.(TPHPOS_BNST_L{a})=struct('BNST',[]);

end

for a=1:length(TPHPOS_BNST_L)
        for b=2%1:length(names)
         forcutoff1_pre.(TPHPOS_BNST_L{a}).(names{b})=cat(2,LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9.meaneachpip_altogether(73:154,:),LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10.meaneachpip_altogether(73:154,:))
    end
end

for a=1:length(TPHPOS_BNST_L)
        for b=2%1:length(names)
        cutoff2_pre.(TPHPOS_BNST_L{a}).(names{b})=median(forcutoff1_pre.(TPHPOS_BNST_L{a}).(names{b}),2)
    end
end

for a=1:length(TPHPOS_BNST_L)
        for b=2%1:length(names)
        cutoff3_pre.(TPHPOS_BNST_L{a}).(names{b})=(mean(cutoff2_pre.(TPHPOS_BNST_L{a}).(names{b})*15));
    end
end

 for a=1:length(TPHPOS_BNST_L)
     for b=2%1:length(names)
          for c=1:length(pCSname)
              for d=1:length(Prepipname)

 	xx_ind1=mean(mean(LFPsall.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d})(73:154,:),2))>(cutoff3_pre.(TPHPOS_BNST_L{a}).(names{b}));
    find(xx_ind1==1)
    if find(xx_ind1==1) ~= 0
     
   LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d})(:,:)=nan;       
    end
              end
          end
     end
 end

%%
% cea
for a=1:length(TPHNEG_CeA_L)
figure (a)
    for b=1 
     subplot(1,2,(b))
        for c=1:10
           
uuuu=LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
plot(mean(uuuu.prepip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHNEG ' char(TPHNEG_CeA_L(a)), 'all pCSs pips ' char(names(b))])


   subplot(1,2,b+1)
        for c=1:10
           
uuuu=LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
plot(mean(uuuu.pip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHNEG ' char(TPHNEG_CeA_L(a)), 'all CSs pips ' char(names(b))])
        end
end


% bnst
for a=1:length(TPHNEG_BNST_L)
figure (a+40)
    for b=2 
     subplot(1,2,(b-1))
        for c=1:10
           
uuuu=LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
plot(mean(uuuu.prepip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHNEG ' char(TPHNEG_BNST_L(a)), 'all pCSs pips ' char(names(b))])


   subplot(1,2,b)
        for c=1:10
           
uuuu=LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
plot(mean(uuuu.pip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHNEG ' char(TPHNEG_BNST_L(a)), 'all CSs pips ' char(names(b))])
        end
end

%% tph+

for a=1:length(TPHPOS_CeA_L)
figure (a)
    for b=1 
     subplot(1,2,(b))
        for c=1:10
           
uuuu=LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
plot(mean(uuuu.prepip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHPOS ' char(TPHPOS_CeA_L(a)), 'all pCSs pips ' char(names(b))])


   subplot(1,2,b+1)
        for c=1:10
           
uuuu=LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
plot(mean(uuuu.pip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHPOS ' char(TPHPOS_CeA_L(a)), 'all CSs pips ' char(names(b))])
        end
end

%bnst
for a=1:length(TPHPOS_BNST_L)
figure (a+50)
    for b=2 
     subplot(1,2,(b-1))
        for c=1:10
           
uuuu=LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
plot(mean(uuuu.prepip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.prepip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHPOS ' char(TPHPOS_BNST_L(a)), 'all pCSs pips ' char(names(b))])

   subplot(1,2,b)
        for c=1:10
           
uuuu=LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
plot(mean(uuuu.pip1,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip2,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip3,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip4,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip5,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip6,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip7,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip8,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip9,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip10,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip11,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip12,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip13,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip14,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip15,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip16,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip17,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip18,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip19,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip20,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip21,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip22,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip23,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip24,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip25,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip26,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip27,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip28,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip29,2)); hold on; xlim([70 150])
plot(mean(uuuu.pip30,2)); hold on; xlim([70 150])
        end
         xlabel('Freq (Hz)')
 title(['TPHPOS ' char(TPHPOS_BNST_L(a)), 'all CSs pips ' char(names(b))])
        end
end

%%
TPHNEG_BNST_L={'TPH645_RECALL','TPH647_RECALL','TPH662_RECALL','TPH891_RECALL','TPH906_RECALL','TPH1006_RECALL','TPH1025_RECALL','TPH1398_RECALL','TPH1406_RECALL'};
TPHPOS_BNST_L={'TPH394_RECALL','TPH655_RECALL','TPH657_RECALL','TPH905_RECALL','TPH907_RECALL','TPH1008_RECALL','TPH1389_RECALL','TPH1390_RECALL','TPH1399_RECALL'};

TPHNEG_CeA_L={'TPH645_RECALL','TPH647_RECALL','TPH891_RECALL','TPH1006_RECALL','TPH1406_RECALL'};
TPHPOS_CeA_L={'TPH655_RECALL','TPH691_RECALL','TPH905_RECALL','TPH907_RECALL','TPH1008_RECALL','TPH1390_RECALL'};

%% 6 Power calculations nan

%TPHNEG
%BNST

for a=1:length(TPHNEG_BNST_L)
    for b=2%1:length(names)
        for c=1:numCS            
            ooooo=LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).altogether=cat(2,ooooo.pip1,ooooo.pip2,ooooo.pip3,ooooo.pip4,ooooo.pip5,ooooo.pip6,ooooo.pip7,ooooo.pip8,ooooo.pip9,ooooo.pip10,ooooo.pip11,ooooo.pip12,ooooo.pip13,ooooo.pip14,ooooo.pip15,ooooo.pip16,ooooo.pip17,ooooo.pip18,ooooo.pip19,ooooo.pip20,ooooo.pip21,ooooo.pip22,ooooo.pip23,ooooo.pip24,ooooo.pip25,ooooo.pip26,ooooo.pip27,ooooo.pip28,ooooo.pip29,ooooo.pip30);
        
            ooooo=LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).altogether=cat(2,ooooo.prepip1,ooooo.prepip2,ooooo.prepip3,ooooo.prepip4,ooooo.prepip5,ooooo.prepip6,ooooo.prepip7,ooooo.prepip8,ooooo.prepip9,ooooo.prepip10,ooooo.prepip11,ooooo.prepip12,ooooo.prepip13,ooooo.prepip14,ooooo.prepip15,ooooo.prepip16,ooooo.prepip17,ooooo.prepip18,ooooo.prepip19,ooooo.prepip20,ooooo.prepip21,ooooo.prepip22,ooooo.prepip23,ooooo.prepip24,ooooo.prepip25,ooooo.prepip26,ooooo.prepip27,ooooo.prepip28,ooooo.prepip29,ooooo.prepip30);
        end
    end
end

for a=1:length(TPHNEG_BNST_L)
    for b=2%1:length(names) 
    ooooo=LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime;
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether,ooooo.CS5.altogether,ooooo.CS6.altogether,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs_mean=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs,2,'omitnan');
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4_mean=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4,2,'omitnan');
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10=cat(2,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10_mean=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10,2,'omitnan');
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2_mean=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10.altogether,2);
    end
end

for a=1:length(TPHNEG_BNST_L)
    for b=2%1:length(names) 
    ooooo=LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime;
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs=cat(2,ooooo.pCS1.altogether,ooooo.pCS2.altogether,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs_mean=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs,2,'omitnan');
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10=cat(2,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10_mean=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10,2,'omitnan');     
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10only=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10.altogether,2);
    end
end

pip_mean_name={'pip1_mean','pip2_mean','pip3_mean','pip4_mean','pip5_mean','pip6_mean','pip7_mean','pip8_mean','pip9_mean','pip10_mean','pip11_mean','pip12_mean','pip13_mean','pip14_mean','pip15_mean','pip16_mean','pip17_mean','pip18_mean','pip19_mean','pip20_mean','pip21_mean','pip22_mean','pip23_mean','pip24_mean','pip25_mean','pip26_mean','pip27_mean','pip28_mean','pip29_mean','pip30_mean'};
prepip_mean_name={'prepip1_mean','prepip2_mean','prepip3_mean','prepip4_mean','prepip5_mean','prepip6_mean','prepip7_mean','prepip8_mean','prepip9_mean','prepip10_mean','prepip11_mean','prepip12_mean','prepip13_mean','prepip14_mean','prepip15_mean','prepip16_mean','prepip17_mean','prepip18_mean','prepip19_mean','prepip20_mean','prepip21_mean','prepip22_mean','prepip23_mean','prepip24_mean','prepip25_mean','prepip26_mean','prepip27_mean','prepip28_mean','prepip29_mean','prepip30_mean'};

for a=1:length(TPHNEG_BNST_L)
    for b=2%1:length(names) 
        for c=1:length(CSname) 
            for d=1:length(Pipname) 
        LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(pip_mean_name{d})=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d}),2,'omitnan');
         end
    end
    end
end


for a=1:length(TPHNEG_BNST_L)
    for b=2%1:length(names) 
        for c=1:length(pCSname) 
            for d=1:length(Prepipname) 
        LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(prepip_mean_name{d})=mean(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d}),2,'omitnan');
         end
    end
    end
end

for a=1:length(TPHNEG_BNST_L)
    for b=2%1:length(names) 
        for c=1:length(CSname) 
            ooooo=LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.pip1_mean,ooooo.pip2_mean,ooooo.pip3_mean,ooooo.pip4_mean,ooooo.pip5_mean,ooooo.pip6_mean,ooooo.pip7_mean,ooooo.pip8_mean,ooooo.pip9_mean,ooooo.pip10_mean,ooooo.pip11_mean,ooooo.pip12_mean,ooooo.pip13_mean,ooooo.pip14_mean,ooooo.pip15_mean,ooooo.pip16_mean,ooooo.pip17_mean,ooooo.pip18_mean,ooooo.pip19_mean,ooooo.pip20_mean,ooooo.pip21_mean,ooooo.pip22_mean,ooooo.pip23_mean,ooooo.pip24_mean,ooooo.pip25_mean,ooooo.pip26_mean,ooooo.pip27_mean,ooooo.pip28_mean,ooooo.pip29_mean,ooooo.pip30_mean);  
        end
    end
end

for a=1:length(TPHNEG_BNST_L)
    for b=2%1:length(names) 
        for c=1:length(pCSname) 
              ooooo=LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.prepip1_mean,ooooo.prepip2_mean,ooooo.prepip3_mean,ooooo.prepip4_mean,ooooo.prepip5_mean,ooooo.prepip6_mean,ooooo.prepip7_mean,ooooo.prepip8_mean,ooooo.prepip9_mean,ooooo.prepip10_mean,ooooo.prepip11_mean,ooooo.prepip12_mean,ooooo.prepip13_mean,ooooo.prepip14_mean,ooooo.prepip15_mean,ooooo.prepip16_mean,ooooo.prepip17_mean,ooooo.prepip18_mean,ooooo.prepip19_mean,ooooo.prepip20_mean,ooooo.prepip21_mean,ooooo.prepip22_mean,ooooo.prepip23_mean,ooooo.prepip24_mean,ooooo.prepip25_mean,ooooo.prepip26_mean,ooooo.prepip27_mean,ooooo.prepip28_mean,ooooo.prepip29_mean,ooooo.prepip30_mean);  
        end
    end
end

%TPHPOS
%BNST

for a=1:length(TPHPOS_BNST_L)
    for b=2%1:length(names)
        for c=1:numCS            
            ooooo=LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).altogether=cat(2,ooooo.pip1,ooooo.pip2,ooooo.pip3,ooooo.pip4,ooooo.pip5,ooooo.pip6,ooooo.pip7,ooooo.pip8,ooooo.pip9,ooooo.pip10,ooooo.pip11,ooooo.pip12,ooooo.pip13,ooooo.pip14,ooooo.pip15,ooooo.pip16,ooooo.pip17,ooooo.pip18,ooooo.pip19,ooooo.pip20,ooooo.pip21,ooooo.pip22,ooooo.pip23,ooooo.pip24,ooooo.pip25,ooooo.pip26,ooooo.pip27,ooooo.pip28,ooooo.pip29,ooooo.pip30);
        
            ooooo=LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).altogether=cat(2,ooooo.prepip1,ooooo.prepip2,ooooo.prepip3,ooooo.prepip4,ooooo.prepip5,ooooo.prepip6,ooooo.prepip7,ooooo.prepip8,ooooo.prepip9,ooooo.prepip10,ooooo.prepip11,ooooo.prepip12,ooooo.prepip13,ooooo.prepip14,ooooo.prepip15,ooooo.prepip16,ooooo.prepip17,ooooo.prepip18,ooooo.prepip19,ooooo.prepip20,ooooo.prepip21,ooooo.prepip22,ooooo.prepip23,ooooo.prepip24,ooooo.prepip25,ooooo.prepip26,ooooo.prepip27,ooooo.prepip28,ooooo.prepip29,ooooo.prepip30);
        end
    end
end

for a=1:length(TPHPOS_BNST_L)
    for b=2%1:length(names) 
    ooooo=LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime;
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether,ooooo.CS5.altogether,ooooo.CS6.altogether,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs_mean=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs,2,'omitnan');
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4_mean=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4,2,'omitnan');
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10=cat(2,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10_mean=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10,2,'omitnan');
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2_mean=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10.altogether,2);
    end
end


for a=1:length(TPHPOS_BNST_L)
    for b=2%1:length(names) 
    ooooo=LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime;
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs=cat(2,ooooo.pCS1.altogether,ooooo.pCS2.altogether,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs_mean=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs,2,'omitnan');
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10=cat(2,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10_mean=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10,2,'omitnan');     
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10only=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10.altogether,2);
    end
end

pip_mean_name={'pip1_mean','pip2_mean','pip3_mean','pip4_mean','pip5_mean','pip6_mean','pip7_mean','pip8_mean','pip9_mean','pip10_mean','pip11_mean','pip12_mean','pip13_mean','pip14_mean','pip15_mean','pip16_mean','pip17_mean','pip18_mean','pip19_mean','pip20_mean','pip21_mean','pip22_mean','pip23_mean','pip24_mean','pip25_mean','pip26_mean','pip27_mean','pip28_mean','pip29_mean','pip30_mean'};
prepip_mean_name={'prepip1_mean','prepip2_mean','prepip3_mean','prepip4_mean','prepip5_mean','prepip6_mean','prepip7_mean','prepip8_mean','prepip9_mean','prepip10_mean','prepip11_mean','prepip12_mean','prepip13_mean','prepip14_mean','prepip15_mean','prepip16_mean','prepip17_mean','prepip18_mean','prepip19_mean','prepip20_mean','prepip21_mean','prepip22_mean','prepip23_mean','prepip24_mean','prepip25_mean','prepip26_mean','prepip27_mean','prepip28_mean','prepip29_mean','prepip30_mean'};

for a=1:length(TPHPOS_BNST_L)
    for b=2%1:length(names) 
        for c=1:length(CSname) 
            for d=1:length(Pipname) 
        LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(pip_mean_name{d})=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d}),2,'omitnan');
         end
    end
    end
end

for a=1:length(TPHPOS_BNST_L)
    for b=2%1:length(names) 
        for c=1:length(pCSname) 
            for d=1:length(Prepipname) 
        LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(prepip_mean_name{d})=mean(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d}),2,'omitnan');
         end
    end
    end
end

for a=1:length(TPHPOS_BNST_L)
    for b=2%1:length(names) 
        for c=1:length(CSname) 
            ooooo=LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.pip1_mean,ooooo.pip2_mean,ooooo.pip3_mean,ooooo.pip4_mean,ooooo.pip5_mean,ooooo.pip6_mean,ooooo.pip7_mean,ooooo.pip8_mean,ooooo.pip9_mean,ooooo.pip10_mean,ooooo.pip11_mean,ooooo.pip12_mean,ooooo.pip13_mean,ooooo.pip14_mean,ooooo.pip15_mean,ooooo.pip16_mean,ooooo.pip17_mean,ooooo.pip18_mean,ooooo.pip19_mean,ooooo.pip20_mean,ooooo.pip21_mean,ooooo.pip22_mean,ooooo.pip23_mean,ooooo.pip24_mean,ooooo.pip25_mean,ooooo.pip26_mean,ooooo.pip27_mean,ooooo.pip28_mean,ooooo.pip29_mean,ooooo.pip30_mean);  
        end
    end
end

for a=1:length(TPHPOS_BNST_L)
    for b=2%1:length(names) 
        for c=1:length(pCSname) 
              ooooo=LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.prepip1_mean,ooooo.prepip2_mean,ooooo.prepip3_mean,ooooo.prepip4_mean,ooooo.prepip5_mean,ooooo.prepip6_mean,ooooo.prepip7_mean,ooooo.prepip8_mean,ooooo.prepip9_mean,ooooo.prepip10_mean,ooooo.prepip11_mean,ooooo.prepip12_mean,ooooo.prepip13_mean,ooooo.prepip14_mean,ooooo.prepip15_mean,ooooo.prepip16_mean,ooooo.prepip17_mean,ooooo.prepip18_mean,ooooo.prepip19_mean,ooooo.prepip20_mean,ooooo.prepip21_mean,ooooo.prepip22_mean,ooooo.prepip23_mean,ooooo.prepip24_mean,ooooo.prepip25_mean,ooooo.prepip26_mean,ooooo.prepip27_mean,ooooo.prepip28_mean,ooooo.prepip29_mean,ooooo.prepip30_mean);  
        end
    end
end

%TPHNEG
%CeA

for a=1:length(TPHNEG_CeA_L)
    for b=1%:length(names)
        for c=1:numCS            
            ooooo=LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).altogether=cat(2,ooooo.pip1,ooooo.pip2,ooooo.pip3,ooooo.pip4,ooooo.pip5,ooooo.pip6,ooooo.pip7,ooooo.pip8,ooooo.pip9,ooooo.pip10,ooooo.pip11,ooooo.pip12,ooooo.pip13,ooooo.pip14,ooooo.pip15,ooooo.pip16,ooooo.pip17,ooooo.pip18,ooooo.pip19,ooooo.pip20,ooooo.pip21,ooooo.pip22,ooooo.pip23,ooooo.pip24,ooooo.pip25,ooooo.pip26,ooooo.pip27,ooooo.pip28,ooooo.pip29,ooooo.pip30);
        
            ooooo=LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).altogether=cat(2,ooooo.prepip1,ooooo.prepip2,ooooo.prepip3,ooooo.prepip4,ooooo.prepip5,ooooo.prepip6,ooooo.prepip7,ooooo.prepip8,ooooo.prepip9,ooooo.prepip10,ooooo.prepip11,ooooo.prepip12,ooooo.prepip13,ooooo.prepip14,ooooo.prepip15,ooooo.prepip16,ooooo.prepip17,ooooo.prepip18,ooooo.prepip19,ooooo.prepip20,ooooo.prepip21,ooooo.prepip22,ooooo.prepip23,ooooo.prepip24,ooooo.prepip25,ooooo.prepip26,ooooo.prepip27,ooooo.prepip28,ooooo.prepip29,ooooo.prepip30);
        end
    end
end

for a=1:length(TPHNEG_CeA_L)
    for b=1%:length(names) 
    ooooo=LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime;
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether,ooooo.CS5.altogether,ooooo.CS6.altogether,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs_mean=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs,2,'omitnan');
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4_mean=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4,2,'omitnan');
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10=cat(2,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10_mean=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10,2,'omitnan');
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2_mean=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10.altogether,2);
    end
end

for a=1:length(TPHNEG_CeA_L)
    for b=1%:length(names) 
    ooooo=LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime;
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs=cat(2,ooooo.pCS1.altogether,ooooo.pCS2.altogether,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs_mean=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs,2,'omitnan');
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10=cat(2,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10_mean=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10,2,'omitnan');     
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9.altogether,2);
    LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10only=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10.altogether,2);
    end
end

pip_mean_name={'pip1_mean','pip2_mean','pip3_mean','pip4_mean','pip5_mean','pip6_mean','pip7_mean','pip8_mean','pip9_mean','pip10_mean','pip11_mean','pip12_mean','pip13_mean','pip14_mean','pip15_mean','pip16_mean','pip17_mean','pip18_mean','pip19_mean','pip20_mean','pip21_mean','pip22_mean','pip23_mean','pip24_mean','pip25_mean','pip26_mean','pip27_mean','pip28_mean','pip29_mean','pip30_mean'};
prepip_mean_name={'prepip1_mean','prepip2_mean','prepip3_mean','prepip4_mean','prepip5_mean','prepip6_mean','prepip7_mean','prepip8_mean','prepip9_mean','prepip10_mean','prepip11_mean','prepip12_mean','prepip13_mean','prepip14_mean','prepip15_mean','prepip16_mean','prepip17_mean','prepip18_mean','prepip19_mean','prepip20_mean','prepip21_mean','prepip22_mean','prepip23_mean','prepip24_mean','prepip25_mean','prepip26_mean','prepip27_mean','prepip28_mean','prepip29_mean','prepip30_mean'};

for a=1:length(TPHNEG_CeA_L)
    for b=1%:length(names) 
        for c=1:length(CSname) 
            for d=1:length(Pipname) 
        LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(pip_mean_name{d})=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d}),2,'omitnan');
         end
    end
    end
end


for a=1:length(TPHNEG_CeA_L)
    for b=1%:length(names) 
        for c=1:length(pCSname) 
            for d=1:length(Prepipname) 
        LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(prepip_mean_name{d})=mean(LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d}),2,'omitnan');
         end
    end
    end
end

for a=1:length(TPHNEG_CeA_L)
    for b=1%:length(names) 
        for c=1:length(CSname) 
            ooooo=LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.pip1_mean,ooooo.pip2_mean,ooooo.pip3_mean,ooooo.pip4_mean,ooooo.pip5_mean,ooooo.pip6_mean,ooooo.pip7_mean,ooooo.pip8_mean,ooooo.pip9_mean,ooooo.pip10_mean,ooooo.pip11_mean,ooooo.pip12_mean,ooooo.pip13_mean,ooooo.pip14_mean,ooooo.pip15_mean,ooooo.pip16_mean,ooooo.pip17_mean,ooooo.pip18_mean,ooooo.pip19_mean,ooooo.pip20_mean,ooooo.pip21_mean,ooooo.pip22_mean,ooooo.pip23_mean,ooooo.pip24_mean,ooooo.pip25_mean,ooooo.pip26_mean,ooooo.pip27_mean,ooooo.pip28_mean,ooooo.pip29_mean,ooooo.pip30_mean);  
        end
    end
end

for a=1:length(TPHNEG_CeA_L)
    for b=1%:length(names) 
        for c=1:length(pCSname) 
              ooooo=LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.prepip1_mean,ooooo.prepip2_mean,ooooo.prepip3_mean,ooooo.prepip4_mean,ooooo.prepip5_mean,ooooo.prepip6_mean,ooooo.prepip7_mean,ooooo.prepip8_mean,ooooo.prepip9_mean,ooooo.prepip10_mean,ooooo.prepip11_mean,ooooo.prepip12_mean,ooooo.prepip13_mean,ooooo.prepip14_mean,ooooo.prepip15_mean,ooooo.prepip16_mean,ooooo.prepip17_mean,ooooo.prepip18_mean,ooooo.prepip19_mean,ooooo.prepip20_mean,ooooo.prepip21_mean,ooooo.prepip22_mean,ooooo.prepip23_mean,ooooo.prepip24_mean,ooooo.prepip25_mean,ooooo.prepip26_mean,ooooo.prepip27_mean,ooooo.prepip28_mean,ooooo.prepip29_mean,ooooo.prepip30_mean);  
        end
    end
end

%TPHPOS
%CeA

for a=1:length(TPHPOS_CeA_L)
    for b=1%:length(names)
        for c=1:numCS            
            ooooo=LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).altogether=cat(2,ooooo.pip1,ooooo.pip2,ooooo.pip3,ooooo.pip4,ooooo.pip5,ooooo.pip6,ooooo.pip7,ooooo.pip8,ooooo.pip9,ooooo.pip10,ooooo.pip11,ooooo.pip12,ooooo.pip13,ooooo.pip14,ooooo.pip15,ooooo.pip16,ooooo.pip17,ooooo.pip18,ooooo.pip19,ooooo.pip20,ooooo.pip21,ooooo.pip22,ooooo.pip23,ooooo.pip24,ooooo.pip25,ooooo.pip26,ooooo.pip27,ooooo.pip28,ooooo.pip29,ooooo.pip30);
        
            ooooo=LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).altogether=cat(2,ooooo.prepip1,ooooo.prepip2,ooooo.prepip3,ooooo.prepip4,ooooo.prepip5,ooooo.prepip6,ooooo.prepip7,ooooo.prepip8,ooooo.prepip9,ooooo.prepip10,ooooo.prepip11,ooooo.prepip12,ooooo.prepip13,ooooo.prepip14,ooooo.prepip15,ooooo.prepip16,ooooo.prepip17,ooooo.prepip18,ooooo.prepip19,ooooo.prepip20,ooooo.prepip21,ooooo.prepip22,ooooo.prepip23,ooooo.prepip24,ooooo.prepip25,ooooo.prepip26,ooooo.prepip27,ooooo.prepip28,ooooo.prepip29,ooooo.prepip30);
        end
    end
end

for a=1:length(TPHPOS_CeA_L)
    for b=1%:length(names) 
    ooooo=LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime;
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether,ooooo.CS5.altogether,ooooo.CS6.altogether,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs_mean=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.all_CSs,2,'omitnan');
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4_mean=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_4,2,'omitnan');
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10=cat(2,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10_mean=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7_10,2,'omitnan');
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2_mean=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1_2,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS1.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS2.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS3.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS4.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS5.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS6.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS7.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS8.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS9.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.CS10.altogether,2);
    end
end

for a=1:length(TPHPOS_CeA_L)
    for b=1%:length(names) 
    ooooo=LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime;
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs=cat(2,ooooo.pCS1.altogether,ooooo.pCS2.altogether,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs_mean=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.all_pCSs,2,'omitnan');
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10=cat(2,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10_mean=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3_10,2,'omitnan');     
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS1.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS2.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS3.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS4.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS5.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS6.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS7.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS8.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS9.altogether,2);
    LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10only=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.pCS10.altogether,2);
    end
end

pip_mean_name={'pip1_mean','pip2_mean','pip3_mean','pip4_mean','pip5_mean','pip6_mean','pip7_mean','pip8_mean','pip9_mean','pip10_mean','pip11_mean','pip12_mean','pip13_mean','pip14_mean','pip15_mean','pip16_mean','pip17_mean','pip18_mean','pip19_mean','pip20_mean','pip21_mean','pip22_mean','pip23_mean','pip24_mean','pip25_mean','pip26_mean','pip27_mean','pip28_mean','pip29_mean','pip30_mean'};
prepip_mean_name={'prepip1_mean','prepip2_mean','prepip3_mean','prepip4_mean','prepip5_mean','prepip6_mean','prepip7_mean','prepip8_mean','prepip9_mean','prepip10_mean','prepip11_mean','prepip12_mean','prepip13_mean','prepip14_mean','prepip15_mean','prepip16_mean','prepip17_mean','prepip18_mean','prepip19_mean','prepip20_mean','prepip21_mean','prepip22_mean','prepip23_mean','prepip24_mean','prepip25_mean','prepip26_mean','prepip27_mean','prepip28_mean','prepip29_mean','prepip30_mean'};

for a=1:length(TPHPOS_CeA_L)
    for b=1%:length(names) 
        for c=1:length(CSname) 
            for d=1:length(Pipname) 
        LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(pip_mean_name{d})=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).(Pipname{d}),2,'omitnan');
         end
    end
    end
end

for a=1:length(TPHPOS_CeA_L)
    for b=1%:length(names) 
        for c=1:length(pCSname) 
            for d=1:length(Prepipname) 
        LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(prepip_mean_name{d})=mean(LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d}),2,'omitnan');
         end
    end
    end
end

for a=1:length(TPHPOS_CeA_L)
    for b=1%:length(names) 
        for c=1:length(CSname) 
            ooooo=LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c});
            LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Pips_freqtime.(CSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.pip1_mean,ooooo.pip2_mean,ooooo.pip3_mean,ooooo.pip4_mean,ooooo.pip5_mean,ooooo.pip6_mean,ooooo.pip7_mean,ooooo.pip8_mean,ooooo.pip9_mean,ooooo.pip10_mean,ooooo.pip11_mean,ooooo.pip12_mean,ooooo.pip13_mean,ooooo.pip14_mean,ooooo.pip15_mean,ooooo.pip16_mean,ooooo.pip17_mean,ooooo.pip18_mean,ooooo.pip19_mean,ooooo.pip20_mean,ooooo.pip21_mean,ooooo.pip22_mean,ooooo.pip23_mean,ooooo.pip24_mean,ooooo.pip25_mean,ooooo.pip26_mean,ooooo.pip27_mean,ooooo.pip28_mean,ooooo.pip29_mean,ooooo.pip30_mean);  
        end
    end
end

for a=1:length(TPHPOS_CeA_L)
    for b=1%:length(names) 
        for c=1:length(pCSname) 
              ooooo=LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c});
            LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{b}).pips.Prepips_freqtime.(pCSname{c}).meaneachpip_altogether(:,:)=cat(2,ooooo.prepip1_mean,ooooo.prepip2_mean,ooooo.prepip3_mean,ooooo.prepip4_mean,ooooo.prepip5_mean,ooooo.prepip6_mean,ooooo.prepip7_mean,ooooo.prepip8_mean,ooooo.prepip9_mean,ooooo.prepip10_mean,ooooo.prepip11_mean,ooooo.prepip12_mean,ooooo.prepip13_mean,ooooo.prepip14_mean,ooooo.prepip15_mean,ooooo.prepip16_mean,ooooo.prepip17_mean,ooooo.prepip18_mean,ooooo.prepip19_mean,ooooo.prepip20_mean,ooooo.prepip21_mean,ooooo.prepip22_mean,ooooo.prepip23_mean,ooooo.prepip24_mean,ooooo.prepip25_mean,ooooo.prepip26_mean,ooooo.prepip27_mean,ooooo.prepip28_mean,ooooo.prepip29_mean,ooooo.prepip30_mean);  
        end
    end
end

 %% pwr calcs

%cea
fooTPHNEG_CeA=LFPsallnan.TPHNEG;
for a=1:length(TPHNEG_CeA_L)
     for n=1%:length(names)
    fooTPHNEG_CeA=LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{n}).pips.Pips_freqtime;
Pwrnan.TPHNEG.(names{n}).CS1_4_mean(:,a)=fooTPHNEG_CeA.CS1_4_mean;
Pwrnan.TPHNEG.(names{n}).CS7_10_mean(:,a)=fooTPHNEG_CeA.CS7_10_mean;
Pwrnan.TPHNEG.(names{n}).all_CSs_mean(:,a)=fooTPHNEG_CeA.all_CSs_mean;
Pwrnan.TPHNEG.(names{n}).CS1_mean(:,a)=fooTPHNEG_CeA.CS1only;
Pwrnan.TPHNEG.(names{n}).CS2_mean(:,a)=fooTPHNEG_CeA.CS2only;
Pwrnan.TPHNEG.(names{n}).CS3_mean(:,a)=fooTPHNEG_CeA.CS3only;
Pwrnan.TPHNEG.(names{n}).CS4_mean(:,a)=fooTPHNEG_CeA.CS4only;
Pwrnan.TPHNEG.(names{n}).CS5_mean(:,a)=fooTPHNEG_CeA.CS5only;
Pwrnan.TPHNEG.(names{n}).CS6_mean(:,a)=fooTPHNEG_CeA.CS6only;
Pwrnan.TPHNEG.(names{n}).CS7_mean(:,a)=fooTPHNEG_CeA.CS7only;
Pwrnan.TPHNEG.(names{n}).CS8_mean(:,a)=fooTPHNEG_CeA.CS8only;
Pwrnan.TPHNEG.(names{n}).CS9_mean(:,a)=fooTPHNEG_CeA.CS9only;
Pwrnan.TPHNEG.(names{n}).CS10_mean(:,a)=fooTPHNEG_CeA.CS10only;

    fooTPHNEG_CeA=LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{n}).pips.Prepips_freqtime; 
Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,a)=fooTPHNEG_CeA.all_pCSs_mean;

Pwrnan.TPHNEG.(names{n}).normCS1_4_mean=Pwrnan.TPHNEG.(names{n}).CS1_4_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).normCS7_10_mean=Pwrnan.TPHNEG.(names{n}).CS7_10_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:); 
Pwrnan.TPHNEG.(names{n}).norm_all_CSs_mean=Pwrnan.TPHNEG.(names{n}).all_CSs_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS1_mean=Pwrnan.TPHNEG.(names{n}).CS1_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS2_mean=Pwrnan.TPHNEG.(names{n}).CS2_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS3_mean=Pwrnan.TPHNEG.(names{n}).CS3_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS4_mean=Pwrnan.TPHNEG.(names{n}).CS4_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS5_mean=Pwrnan.TPHNEG.(names{n}).CS5_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS6_mean=Pwrnan.TPHNEG.(names{n}).CS6_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS7_mean=Pwrnan.TPHNEG.(names{n}).CS7_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS8_mean=Pwrnan.TPHNEG.(names{n}).CS8_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS9_mean=Pwrnan.TPHNEG.(names{n}).CS9_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS10_mean=Pwrnan.TPHNEG.(names{n}).CS10_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);

end
end

fooTPHPOS_CeA=LFPsallnan.TPHPOS;
for a=1:length(TPHPOS_CeA_L)
     for n=1%:length(names)
    fooTPHPOS_CeA=LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{n}).pips.Pips_freqtime;
Pwrnan.TPHPOS.(names{n}).CS1_4_mean(:,a)=fooTPHPOS_CeA.CS1_4_mean;
Pwrnan.TPHPOS.(names{n}).CS7_10_mean(:,a)=fooTPHPOS_CeA.CS7_10_mean;
Pwrnan.TPHPOS.(names{n}).all_CSs_mean(:,a)=fooTPHPOS_CeA.all_CSs_mean;
Pwrnan.TPHPOS.(names{n}).CS1_mean(:,a)=fooTPHPOS_CeA.CS1only;
Pwrnan.TPHPOS.(names{n}).CS2_mean(:,a)=fooTPHPOS_CeA.CS2only;
Pwrnan.TPHPOS.(names{n}).CS3_mean(:,a)=fooTPHPOS_CeA.CS3only;
Pwrnan.TPHPOS.(names{n}).CS4_mean(:,a)=fooTPHPOS_CeA.CS4only;
Pwrnan.TPHPOS.(names{n}).CS5_mean(:,a)=fooTPHPOS_CeA.CS5only;
Pwrnan.TPHPOS.(names{n}).CS6_mean(:,a)=fooTPHPOS_CeA.CS6only;
Pwrnan.TPHPOS.(names{n}).CS7_mean(:,a)=fooTPHPOS_CeA.CS7only;
Pwrnan.TPHPOS.(names{n}).CS8_mean(:,a)=fooTPHPOS_CeA.CS8only;
Pwrnan.TPHPOS.(names{n}).CS9_mean(:,a)=fooTPHPOS_CeA.CS9only;
Pwrnan.TPHPOS.(names{n}).CS10_mean(:,a)=fooTPHPOS_CeA.CS10only;

    fooTPHPOS_CeA=LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{n}).pips.Prepips_freqtime; 
Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,a)=fooTPHPOS_CeA.all_pCSs_mean;

Pwrnan.TPHPOS.(names{n}).normCS1_4_mean=Pwrnan.TPHPOS.(names{n}).CS1_4_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).normCS7_10_mean=Pwrnan.TPHPOS.(names{n}).CS7_10_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:); 
Pwrnan.TPHPOS.(names{n}).norm_all_CSs_mean=Pwrnan.TPHPOS.(names{n}).all_CSs_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS1_mean=Pwrnan.TPHPOS.(names{n}).CS1_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS2_mean=Pwrnan.TPHPOS.(names{n}).CS2_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS3_mean=Pwrnan.TPHPOS.(names{n}).CS3_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS4_mean=Pwrnan.TPHPOS.(names{n}).CS4_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS5_mean=Pwrnan.TPHPOS.(names{n}).CS5_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS6_mean=Pwrnan.TPHPOS.(names{n}).CS6_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS7_mean=Pwrnan.TPHPOS.(names{n}).CS7_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS8_mean=Pwrnan.TPHPOS.(names{n}).CS8_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS9_mean=Pwrnan.TPHPOS.(names{n}).CS9_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS10_mean=Pwrnan.TPHPOS.(names{n}).CS10_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);

end
end

%bnst
fooTPHNEG_BNST=LFPsallnan.TPHNEG;
for a=1:length(TPHNEG_BNST_L)
     for n=2%1:length(names)
    fooTPHNEG_BNST=LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{n}).pips.Pips_freqtime;
Pwrnan.TPHNEG.(names{n}).CS1_4_mean(:,a)=fooTPHNEG_BNST.CS1_4_mean;
Pwrnan.TPHNEG.(names{n}).CS7_10_mean(:,a)=fooTPHNEG_BNST.CS7_10_mean;
Pwrnan.TPHNEG.(names{n}).all_CSs_mean(:,a)=fooTPHNEG_BNST.all_CSs_mean;
Pwrnan.TPHNEG.(names{n}).CS1_mean(:,a)=fooTPHNEG_BNST.CS1only;
Pwrnan.TPHNEG.(names{n}).CS2_mean(:,a)=fooTPHNEG_BNST.CS2only;
Pwrnan.TPHNEG.(names{n}).CS3_mean(:,a)=fooTPHNEG_BNST.CS3only;
Pwrnan.TPHNEG.(names{n}).CS4_mean(:,a)=fooTPHNEG_BNST.CS4only;
Pwrnan.TPHNEG.(names{n}).CS5_mean(:,a)=fooTPHNEG_BNST.CS5only;
Pwrnan.TPHNEG.(names{n}).CS6_mean(:,a)=fooTPHNEG_BNST.CS6only;
Pwrnan.TPHNEG.(names{n}).CS7_mean(:,a)=fooTPHNEG_BNST.CS7only;
Pwrnan.TPHNEG.(names{n}).CS8_mean(:,a)=fooTPHNEG_BNST.CS8only;
Pwrnan.TPHNEG.(names{n}).CS9_mean(:,a)=fooTPHNEG_BNST.CS9only;
Pwrnan.TPHNEG.(names{n}).CS10_mean(:,a)=fooTPHNEG_BNST.CS10only;

    fooTPHNEG_BNST=LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{n}).pips.Prepips_freqtime; 
Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,a)=fooTPHNEG_BNST.all_pCSs_mean;

Pwrnan.TPHNEG.(names{n}).normCS1_4_mean=Pwrnan.TPHNEG.(names{n}).CS1_4_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).normCS7_10_mean=Pwrnan.TPHNEG.(names{n}).CS7_10_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:); 
Pwrnan.TPHNEG.(names{n}).norm_all_CSs_mean=Pwrnan.TPHNEG.(names{n}).all_CSs_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS1_mean=Pwrnan.TPHNEG.(names{n}).CS1_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS2_mean=Pwrnan.TPHNEG.(names{n}).CS2_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS3_mean=Pwrnan.TPHNEG.(names{n}).CS3_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS4_mean=Pwrnan.TPHNEG.(names{n}).CS4_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS5_mean=Pwrnan.TPHNEG.(names{n}).CS5_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS6_mean=Pwrnan.TPHNEG.(names{n}).CS6_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS7_mean=Pwrnan.TPHNEG.(names{n}).CS7_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS8_mean=Pwrnan.TPHNEG.(names{n}).CS8_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS9_mean=Pwrnan.TPHNEG.(names{n}).CS9_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS10_mean=Pwrnan.TPHNEG.(names{n}).CS10_mean(:,:)./Pwrnan.TPHNEG.(names{n}).all_pCSs_mean(:,:);

end
end

fooTPHPOS_BNST=LFPsallnan.TPHPOS;
for a=1:length(TPHPOS_BNST_L)
     for n=2%1:length(names)
    fooTPHPOS_BNST=LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{n}).pips.Pips_freqtime;
Pwrnan.TPHPOS.(names{n}).CS1_4_mean(:,a)=fooTPHPOS_BNST.CS1_4_mean;
Pwrnan.TPHPOS.(names{n}).CS7_10_mean(:,a)=fooTPHPOS_BNST.CS7_10_mean;
Pwrnan.TPHPOS.(names{n}).all_CSs_mean(:,a)=fooTPHPOS_BNST.all_CSs_mean;
Pwrnan.TPHPOS.(names{n}).CS1_mean(:,a)=fooTPHPOS_BNST.CS1only;
Pwrnan.TPHPOS.(names{n}).CS2_mean(:,a)=fooTPHPOS_BNST.CS2only;
Pwrnan.TPHPOS.(names{n}).CS3_mean(:,a)=fooTPHPOS_BNST.CS3only;
Pwrnan.TPHPOS.(names{n}).CS4_mean(:,a)=fooTPHPOS_BNST.CS4only;
Pwrnan.TPHPOS.(names{n}).CS5_mean(:,a)=fooTPHPOS_BNST.CS5only;
Pwrnan.TPHPOS.(names{n}).CS6_mean(:,a)=fooTPHPOS_BNST.CS6only;
Pwrnan.TPHPOS.(names{n}).CS7_mean(:,a)=fooTPHPOS_BNST.CS7only;
Pwrnan.TPHPOS.(names{n}).CS8_mean(:,a)=fooTPHPOS_BNST.CS8only;
Pwrnan.TPHPOS.(names{n}).CS9_mean(:,a)=fooTPHPOS_BNST.CS9only;
Pwrnan.TPHPOS.(names{n}).CS10_mean(:,a)=fooTPHPOS_BNST.CS10only;

    fooTPHPOS_BNST=LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{n}).pips.Prepips_freqtime; 
Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,a)=fooTPHPOS_BNST.all_pCSs_mean;

Pwrnan.TPHPOS.(names{n}).normCS1_4_mean=Pwrnan.TPHPOS.(names{n}).CS1_4_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).normCS7_10_mean=Pwrnan.TPHPOS.(names{n}).CS7_10_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:); 
Pwrnan.TPHPOS.(names{n}).norm_all_CSs_mean=Pwrnan.TPHPOS.(names{n}).all_CSs_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS1_mean=Pwrnan.TPHPOS.(names{n}).CS1_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS2_mean=Pwrnan.TPHPOS.(names{n}).CS2_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS3_mean=Pwrnan.TPHPOS.(names{n}).CS3_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS4_mean=Pwrnan.TPHPOS.(names{n}).CS4_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS5_mean=Pwrnan.TPHPOS.(names{n}).CS5_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS6_mean=Pwrnan.TPHPOS.(names{n}).CS6_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS7_mean=Pwrnan.TPHPOS.(names{n}).CS7_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS8_mean=Pwrnan.TPHPOS.(names{n}).CS8_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS9_mean=Pwrnan.TPHPOS.(names{n}).CS9_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS10_mean=Pwrnan.TPHPOS.(names{n}).CS10_mean(:,:)./Pwrnan.TPHPOS.(names{n}).all_pCSs_mean(:,:);

end
end

%%
%Absolute Power for each region, CS and pCS
what={'all_CSs_mean', 'all_pCSs_mean','CS1_4_mean','CS7_10_mean'};

for n=1:length(names)
    for t=1:length(what)
     Pwrnan_all.(names{n}).TPHPOS.(what{t})=(Pwrnan.TPHPOS.(names{n}).(what{t}));
     Pwrnan_all.(names{n}).TPHNEG.(what{t})=(Pwrnan.TPHNEG.(names{n}).(what{t}));
    end
end

Freq=Freq_pw;
gamma_freq_idx_90_140=find(Freq>=90 & Freq<=140);

figure
n=1 % cea

 subplot(3,4,1)     
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{1})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]); ylim([1 100000]) ;title([char(names(n)), ' RECALL TPH-' ]); 
        
        subplot(3,4,5)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{3})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 100000]); title([char(names(n)), ' RECALL TPH-' ]);
        
        subplot(3,4,9)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{4})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 100000]); title([char(names(n)), ' RECALL TPH-']); 

%TPH+ ACSF
     subplot(3,4,2)     
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{1})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]); ylim([1 100000]) ;title([char(names(n)), ' RECALL TPH+' ]); 
        
        subplot(3,4,6)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{3})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0   0.5000])
        xlim([70 150]);ylim([1 100000]); title([char(names(n)), ' RECALL TPH+' ]);
        
        subplot(3,4,10)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{4})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 100000]); title([char(names(n)), ' RECALL TPH+']); 
       
        %TPH- ACSF
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        %CeA L
        subplot(3,4,3)
        CS_SEM = std(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH- all CSs ', ' ','TPH- all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,7)
        CS_SEM = std(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
         ylim([0 80000]);legend('TPH- CS1-4 ', ' ','TPH- all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,11)
        CS_SEM = std(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH- CS7-10 ', ' ','TPH- all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
   
        %TPH+ acsf
       fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        %CeA L
        subplot(3,4,4)
        CS_SEM = std(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH+ all CSs ', ' ','TPH+ all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,8)
        CS_SEM = std(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
         ylim([0 80000]);legend('TPH+ CS1-4 ', ' ','TPH+ all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,12)
        CS_SEM = std(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH+ CS7-10 ', ' ','TPH+ all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])

%%
figure
n=2 % BNST

 subplot(3,4,1)     
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{1})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]); ylim([1 60000]) ;title([char(names(n)), ' RECALL TPH-' ]); 
        
        subplot(3,4,5)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{3})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 60000]); title([char(names(n)), ' RECALL TPH-' ]);
        
        subplot(3,4,9)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{4})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 60000]); title([char(names(n)), ' RECALL TPH-']); 

%TPH+ ACSF
     subplot(3,4,2)     
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{1})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]); ylim([1 60000]) ;title([char(names(n)), ' RECALL TPH+' ]); 
        
        subplot(3,4,6)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{3})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 60000]); title([char(names(n)), ' RECALL TPH+' ]);
        
        subplot(3,4,10)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{4})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 60000]); title([char(names(n)), ' RECALL TPH+']); 
       
        %TPH- ACSF
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        %CeA L
        subplot(3,4,3)
        CS_SEM = std(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 50000]);legend('TPH- all CSs ', ' ','TPH- all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,7)
        CS_SEM = std(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
         ylim([0 50000]);legend('TPH- CS1-4 ', ' ','TPH- all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,11)
        CS_SEM = std(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 50000]);legend('TPH- CS7-10 ', ' ','TPH- all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
   
        %TPH+ acsf
       fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        %CeA L
        subplot(3,4,4)
        CS_SEM = std(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 50000]);legend('TPH+ all CSs ', ' ','TPH+ all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,8)
        CS_SEM = std(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
         ylim([0 50000]);legend('TPH+ CS1-4 ', ' ','TPH+ all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,12)
        CS_SEM = std(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 50000]);legend('TPH+ CS7-10 ', ' ','TPH+ all pCSs'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])

%%
%NORM Power for each region, CS and pCS
clear what
clear Pwrnan_all
what={'norm_all_CSs_mean', 'normCS1_4_mean', 'normCS7_10_mean'};

for n=1:length(names)
    for t=1:length(what)
     Pwrnan_all.(names{n}).TPHPOS.(what{t})=(Pwrnan.TPHPOS.(names{n}).(what{t}));
     Pwrnan_all.(names{n}).TPHNEG.(what{t})=(Pwrnan.TPHNEG.(names{n}).(what{t}));
    end
end

figure(23)   
n=1 ;%cea
        subplot(3,4,1)     
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{1})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{1})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')
       
        subplot(3,4,5)  
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')
        
        subplot(3,4,9)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{3})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{3})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')

subplot(3,4,2)
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')  

        fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1.5]);legend('TPH- all CS ',' ','TPH+ all CS ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
      
        subplot(3,4,6)
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
  
        fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
           ylim([0 1.5]);legend('TPH- CS1-4 ',' ','TPH+ CS1-4 ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
      
        subplot(3,4,10)
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        

        fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1.5]);legend('TPH- CS7-10 ', ' ','TPH+ CS7-10 ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])      

n=2 ;%bnst
        subplot(3,4,3)     
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{1})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{1})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')
       
        subplot(3,4,7)  
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')
        
        subplot(3,4,11)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{3})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{3})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')

subplot(3,4,4)
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')  

        fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1.5]);legend('TPH- all CS ',' ','TPH+ all CS ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
      
        subplot(3,4,8)
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
  
        fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
           ylim([0 1.5]);legend('TPH- CS1-4 ',' ','TPH+ CS1-4 ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
      
        subplot(3,4,12)
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        
        fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1.5]);legend('TPH- CS7-10 ', ' ','TPH+ CS7-10 ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])      

%%

[cwt,t,f]=cwt(Spectrograms_groups.TPHPOS.(whatp{v}).CS1_4_norm_mean,1:0.1:12,3,2000);

%diff=LFPsallnan.TPHNEG.TPH645_RECALL.unFilt.BNST_L.CS{1, 1}(2001:32000)-LFPsallnan.TPHNEG.TPH645_RECALL.unFilt.BNST_L.CS{1, 10}(2001:32000);

for i=1:length(TPHNEG_BNST_L)
clim=[0 30];
figure;
subplot(3,1,1)
[cw1, t, f]=cwt(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{i}).unFilt.BNST_L.CS{1, 1}(2001:32000),15:.1:150,3,2000);
power1=abs(cw1);
imagesc(t,f,power1,clim); axis xy; title([char(TPHPOS_BNST_L(i)) ' BNST CS1 unfil '])
colorbar

subplot(3,1,2)
[cw10, t, f]=cwt(LFPsallnan.TPHNEG.(TPHNEG_BNST_L{i}).unFilt.BNST_L.CS{1, 10}(2001:32000),15:.1:150,3,2000);
power10=abs(cw10);
imagesc(t,f,power10,clim); axis xy; title([char(TPHPOS_BNST_L(i)) ' BNST CS10 unfil '])
colorbar

subplot(3,1,3)
clim=[0 15];
%[cw3, t, f]=cwt(diff,15:.1:150,3,2000);
power10_1=power10-power1;
imagesc(t,f,power10_1,clim); axis xy; title([char(TPHPOS_BNST_L(i)) ' BNST CS10-CS1 unfil '])
colorbar

end

for i=1:length(TPHPOS_BNST_L)
clim=[0 30];
figure;
subplot(3,1,1)
[cw1, t, f]=cwt(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{i}).unFilt.BNST_L.CS{1, 1}(2001:32000),15:.1:150,3,2000);
power1=abs(cw1);
imagesc(t,f,power1,clim); axis xy; title([char(TPHPOS_BNST_L(i)) ' BNST+ CS1 unfil '])
colorbar

subplot(3,1,2)
[cw10, t, f]=cwt(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{i}).unFilt.BNST_L.CS{1, 10}(2001:32000),15:.1:150,3,2000);
power10=abs(cw10);
imagesc(t,f,power10,clim); axis xy; title([char(TPHPOS_BNST_L(i)) '+ BNST CS10 unfil '])
colorbar

subplot(3,1,3)
clim=[0 15];
%[cw3, t, f]=cwt(diff,15:.1:150,3,2000);
power10_1=power10-power1;
imagesc(t,f,power10_1,clim); axis xy; title([char(TPHPOS_BNST_L(i)) '+ BNST CS10-CS1 unfil '])
colorbar

end


ColorPlot(Time_pw(5:155),Freq_pw,Spectrograms_groups.TPHPOS.(whatp{v}).CS1_4_norm_mean)
ylim([0 12]); xlim([1 15.8]) ;
 xlabel('Time within tone')
 colorbar; 
 ylabel('Frequency (Hz)')
 title([ 'FEMALES Recall BNST Tone-Bin 1-2 (normalized)']) 
 xlabel(' Time since tone onset (s)')
xticks([0 2 4 6 8 10 12 14]);
xticklabels({'0' '2' '4' '6' '8' '10' '12' '14'})
  set(gcf,'color','white')
        set(gca,'FontSize',11)
        set(gca,'fontname','arial')


        %%
for i=5%1:length(TPHPOS_BNST_L)
clim=[0 30];
figure;
subplot(3,1,1)
[cw1, t, f]=cwt(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{i}).unFilt.BNST_L.CS{1, 1}(2001:32000),1:.1:150,3,2000);
power1=abs(cw1);
imagesc(t,f,power1,clim); axis xy; title([char(TPHPOS_BNST_L(i)) ' BNST+ CS1 unfil ']); ylim([1 40]);
colorbar

subplot(3,1,2)
[cw10, t, f]=cwt(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{i}).unFilt.BNST_L.CS{1, 10}(2001:32000),1:.1:150,3,2000);
power10=abs(cw10);
imagesc(t,f,power10,clim); axis xy; title([char(TPHPOS_BNST_L(i)) '+ BNST CS10 unfil ']); ylim([1 40]);
colorbar

subplot(3,1,3)
clim=[0 15];
%[cw3, t, f]=cwt(diff,15:.1:150,3,2000);
power10_1=power10-power1;
imagesc(t,f,power10_1,clim); axis xy; title([char(TPHPOS_BNST_L(i)) '+ BNST CS10-CS1 unfil ']); ylim([1 40]);
colorbar

end

        %%
for i=5%1:length(TPHPOS_BNST_L)
clim=[0 15];
figure;
subplot(3,1,1)
[cw1, t, f]=cwt(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{i}).unFilt.BNST_L.CS{1, 1}(2001:32000),1:.1:150,3,2000);
power1=abs(cw1);
imagesc(t,f,power1,clim); axis xy; title([char(TPHPOS_BNST_L(i)) ' BNST+ CS1 unfil ']); ylim([70 150]);
colorbar

subplot(3,1,2)
[cw10, t, f]=cwt(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{i}).unFilt.BNST_L.CS{1, 10}(2001:32000),1:.1:150,3,2000);
power10=abs(cw10);
imagesc(t,f,power10,clim); axis xy; title([char(TPHPOS_BNST_L(i)) '+ BNST CS10 unfil ']); ylim([70 150]);
colorbar

subplot(3,1,3)
clim=[0 8];
%[cw3, t, f]=cwt(diff,15:.1:150,3,2000);
power10_1=power10-power1;
imagesc(t,f,power10_1,clim); axis xy; title([char(TPHPOS_BNST_L(i)) '+ BNST CS10-CS1 unfil ']); ylim([70 150]);
colorbar
end

%% norm

for i=5%1:length(TPHPOS_BNST_L)
clim=[0 5];
figure;
subplot(3,1,1)
[cwCS1, t, f]=cwt(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{i}).unFilt.BNST_L.CS{1, 1}(2001:32000),1:.1:150,3,2000);
powerCS1=abs(cwCS1);
[cwpCS1, t, f]=cwt(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{i}).unFilt.BNST_L.pCS{1, 1}(2001:32000),1:.1:150,3,2000);
powerpCS1=abs(cwpCS1);
normpower1=powerCS1./powerpCS1;
imagesc(t,f,normpower1,clim); axis xy; title([char(TPHPOS_BNST_L(i)) ' BNST+ normCS1 unfil ']); ylim([12 20]);
colorbar

subplot(3,1,2)
[cwCS10, t, f]=cwt(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{i}).unFilt.BNST_L.CS{1, 10}(2001:32000),1:.1:150,3,2000);
powerCS10=abs(cwCS10);
normpower10=powerCS10./powerpCS1;
imagesc(t,f,normpower10,clim); axis xy; title([char(TPHPOS_BNST_L(i)) ' BNST+ normCS10 unfil ']); ylim([12 20]);
colorbar

subplot(3,1,3)
clim=[0 1];
%[cw3, t, f]=cwt(diff,15:.1:150,3,2000);
normpower10_1=normpower10-normpower1;
imagesc(t,f,normpower10_1,clim); axis xy; title([char(TPHPOS_BNST_L(i)) '+ BNST normCS10-CS1 unfil ']); ylim([12 20]);
colorbar
end

for i=5%1:length(TPHPOS_BNST_L)
clim=[0 3];
figure;
subplot(3,1,1)
[cwCS1, t, f]=cwt(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{i}).unFilt.BNST_L.CS{1, 1}(2001:32000),1:.1:150,3,2000);
powerCS1=abs(cwCS1);
[cwpCS1, t, f]=cwt(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{i}).unFilt.BNST_L.pCS{1, 1}(2001:32000),1:.1:150,3,2000);
powerpCS1=abs(cwpCS1);
normpower1=powerCS1./powerpCS1;
imagesc(t,f,normpower1,clim); axis xy; title([char(TPHPOS_BNST_L(i)) ' BNST+ normCS1 unfil ']); ylim([20 50]);
colorbar

subplot(3,1,2)
[cwCS10, t, f]=cwt(LFPsallnan.TPHPOS.(TPHPOS_BNST_L{i}).unFilt.BNST_L.CS{1, 10}(2001:32000),1:.1:150,3,2000);
powerCS10=abs(cwCS10);
normpower10=powerCS10./powerpCS1;
imagesc(t,f,normpower10,clim); axis xy; title([char(TPHPOS_BNST_L(i)) ' BNST+ normCS10 unfil ']); ylim([20 50]);
colorbar

subplot(3,1,3)
clim=[0 1];
%[cw3, t, f]=cwt(diff,15:.1:150,3,2000);
normpower10_1=normpower10-normpower1;
imagesc(t,f,normpower10_1,clim); axis xy; title([char(TPHPOS_BNST_L(i)) '+ BNST normCS10-CS1 unfil ']); ylim([20 50]);
colorbar
end

%% now if we normalize by CS3-10
%cea
fooTPHNEG_CeA=LFPsallnan.TPHNEG;
for a=1:length(TPHNEG_CeA_L)
     for n=1%:length(names)
    fooTPHNEG_CeA=LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{n}).pips.Pips_freqtime;
Pwrnan.TPHNEG.(names{n}).CS1_4_mean(:,a)=fooTPHNEG_CeA.CS1_4_mean;
Pwrnan.TPHNEG.(names{n}).CS7_10_mean(:,a)=fooTPHNEG_CeA.CS7_10_mean;
Pwrnan.TPHNEG.(names{n}).all_CSs_mean(:,a)=fooTPHNEG_CeA.all_CSs_mean;
Pwrnan.TPHNEG.(names{n}).CS1_mean(:,a)=fooTPHNEG_CeA.CS1only;
Pwrnan.TPHNEG.(names{n}).CS2_mean(:,a)=fooTPHNEG_CeA.CS2only;
Pwrnan.TPHNEG.(names{n}).CS3_mean(:,a)=fooTPHNEG_CeA.CS3only;
Pwrnan.TPHNEG.(names{n}).CS4_mean(:,a)=fooTPHNEG_CeA.CS4only;
Pwrnan.TPHNEG.(names{n}).CS5_mean(:,a)=fooTPHNEG_CeA.CS5only;
Pwrnan.TPHNEG.(names{n}).CS6_mean(:,a)=fooTPHNEG_CeA.CS6only;
Pwrnan.TPHNEG.(names{n}).CS7_mean(:,a)=fooTPHNEG_CeA.CS7only;
Pwrnan.TPHNEG.(names{n}).CS8_mean(:,a)=fooTPHNEG_CeA.CS8only;
Pwrnan.TPHNEG.(names{n}).CS9_mean(:,a)=fooTPHNEG_CeA.CS9only;
Pwrnan.TPHNEG.(names{n}).CS10_mean(:,a)=fooTPHNEG_CeA.CS10only;

    fooTPHNEG_CeA=LFPsallnan.TPHNEG.(TPHNEG_CeA_L{a}).TonePwr.(names{n}).pips.Prepips_freqtime; 
Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,a)=fooTPHNEG_CeA.pCS3_10_mean;

Pwrnan.TPHNEG.(names{n}).normCS1_4_mean=Pwrnan.TPHNEG.(names{n}).CS1_4_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).normCS7_10_mean=Pwrnan.TPHNEG.(names{n}).CS7_10_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:); 
Pwrnan.TPHNEG.(names{n}).norm_all_CSs_mean=Pwrnan.TPHNEG.(names{n}).all_CSs_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS1_mean=Pwrnan.TPHNEG.(names{n}).CS1_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS2_mean=Pwrnan.TPHNEG.(names{n}).CS2_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS3_mean=Pwrnan.TPHNEG.(names{n}).CS3_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS4_mean=Pwrnan.TPHNEG.(names{n}).CS4_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS5_mean=Pwrnan.TPHNEG.(names{n}).CS5_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS6_mean=Pwrnan.TPHNEG.(names{n}).CS6_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS7_mean=Pwrnan.TPHNEG.(names{n}).CS7_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS8_mean=Pwrnan.TPHNEG.(names{n}).CS8_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS9_mean=Pwrnan.TPHNEG.(names{n}).CS9_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS10_mean=Pwrnan.TPHNEG.(names{n}).CS10_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);

end
end

fooTPHPOS_CeA=LFPsallnan.TPHPOS;
for a=1:length(TPHPOS_CeA_L)
     for n=1%:length(names)
    fooTPHPOS_CeA=LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{n}).pips.Pips_freqtime;
Pwrnan.TPHPOS.(names{n}).CS1_4_mean(:,a)=fooTPHPOS_CeA.CS1_4_mean;
Pwrnan.TPHPOS.(names{n}).CS7_10_mean(:,a)=fooTPHPOS_CeA.CS7_10_mean;
Pwrnan.TPHPOS.(names{n}).all_CSs_mean(:,a)=fooTPHPOS_CeA.all_CSs_mean;
Pwrnan.TPHPOS.(names{n}).CS1_mean(:,a)=fooTPHPOS_CeA.CS1only;
Pwrnan.TPHPOS.(names{n}).CS2_mean(:,a)=fooTPHPOS_CeA.CS2only;
Pwrnan.TPHPOS.(names{n}).CS3_mean(:,a)=fooTPHPOS_CeA.CS3only;
Pwrnan.TPHPOS.(names{n}).CS4_mean(:,a)=fooTPHPOS_CeA.CS4only;
Pwrnan.TPHPOS.(names{n}).CS5_mean(:,a)=fooTPHPOS_CeA.CS5only;
Pwrnan.TPHPOS.(names{n}).CS6_mean(:,a)=fooTPHPOS_CeA.CS6only;
Pwrnan.TPHPOS.(names{n}).CS7_mean(:,a)=fooTPHPOS_CeA.CS7only;
Pwrnan.TPHPOS.(names{n}).CS8_mean(:,a)=fooTPHPOS_CeA.CS8only;
Pwrnan.TPHPOS.(names{n}).CS9_mean(:,a)=fooTPHPOS_CeA.CS9only;
Pwrnan.TPHPOS.(names{n}).CS10_mean(:,a)=fooTPHPOS_CeA.CS10only;

    fooTPHPOS_CeA=LFPsallnan.TPHPOS.(TPHPOS_CeA_L{a}).TonePwr.(names{n}).pips.Prepips_freqtime; 
Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,a)=fooTPHPOS_CeA.pCS3_10_mean;

Pwrnan.TPHPOS.(names{n}).normCS1_4_mean=Pwrnan.TPHPOS.(names{n}).CS1_4_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).normCS7_10_mean=Pwrnan.TPHPOS.(names{n}).CS7_10_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:); 
Pwrnan.TPHPOS.(names{n}).norm_all_CSs_mean=Pwrnan.TPHPOS.(names{n}).all_CSs_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS1_mean=Pwrnan.TPHPOS.(names{n}).CS1_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS2_mean=Pwrnan.TPHPOS.(names{n}).CS2_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS3_mean=Pwrnan.TPHPOS.(names{n}).CS3_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS4_mean=Pwrnan.TPHPOS.(names{n}).CS4_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS5_mean=Pwrnan.TPHPOS.(names{n}).CS5_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS6_mean=Pwrnan.TPHPOS.(names{n}).CS6_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS7_mean=Pwrnan.TPHPOS.(names{n}).CS7_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS8_mean=Pwrnan.TPHPOS.(names{n}).CS8_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS9_mean=Pwrnan.TPHPOS.(names{n}).CS9_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS10_mean=Pwrnan.TPHPOS.(names{n}).CS10_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);

end
end

%bnst
fooTPHNEG_BNST=LFPsallnan.TPHNEG;
for a=1:length(TPHNEG_BNST_L)
     for n=2%1:length(names)
    fooTPHNEG_BNST=LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{n}).pips.Pips_freqtime;
Pwrnan.TPHNEG.(names{n}).CS1_4_mean(:,a)=fooTPHNEG_BNST.CS1_4_mean;
Pwrnan.TPHNEG.(names{n}).CS7_10_mean(:,a)=fooTPHNEG_BNST.CS7_10_mean;
Pwrnan.TPHNEG.(names{n}).all_CSs_mean(:,a)=fooTPHNEG_BNST.all_CSs_mean;
Pwrnan.TPHNEG.(names{n}).CS1_mean(:,a)=fooTPHNEG_BNST.CS1only;
Pwrnan.TPHNEG.(names{n}).CS2_mean(:,a)=fooTPHNEG_BNST.CS2only;
Pwrnan.TPHNEG.(names{n}).CS3_mean(:,a)=fooTPHNEG_BNST.CS3only;
Pwrnan.TPHNEG.(names{n}).CS4_mean(:,a)=fooTPHNEG_BNST.CS4only;
Pwrnan.TPHNEG.(names{n}).CS5_mean(:,a)=fooTPHNEG_BNST.CS5only;
Pwrnan.TPHNEG.(names{n}).CS6_mean(:,a)=fooTPHNEG_BNST.CS6only;
Pwrnan.TPHNEG.(names{n}).CS7_mean(:,a)=fooTPHNEG_BNST.CS7only;
Pwrnan.TPHNEG.(names{n}).CS8_mean(:,a)=fooTPHNEG_BNST.CS8only;
Pwrnan.TPHNEG.(names{n}).CS9_mean(:,a)=fooTPHNEG_BNST.CS9only;
Pwrnan.TPHNEG.(names{n}).CS10_mean(:,a)=fooTPHNEG_BNST.CS10only;

    fooTPHNEG_BNST=LFPsallnan.TPHNEG.(TPHNEG_BNST_L{a}).TonePwr.(names{n}).pips.Prepips_freqtime; 
Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,a)=fooTPHNEG_BNST.pCS3_10_mean;

Pwrnan.TPHNEG.(names{n}).normCS1_4_mean=Pwrnan.TPHNEG.(names{n}).CS1_4_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).normCS7_10_mean=Pwrnan.TPHNEG.(names{n}).CS7_10_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:); 
Pwrnan.TPHNEG.(names{n}).norm_all_CSs_mean=Pwrnan.TPHNEG.(names{n}).all_CSs_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS1_mean=Pwrnan.TPHNEG.(names{n}).CS1_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS2_mean=Pwrnan.TPHNEG.(names{n}).CS2_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS3_mean=Pwrnan.TPHNEG.(names{n}).CS3_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS4_mean=Pwrnan.TPHNEG.(names{n}).CS4_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS5_mean=Pwrnan.TPHNEG.(names{n}).CS5_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS6_mean=Pwrnan.TPHNEG.(names{n}).CS6_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS7_mean=Pwrnan.TPHNEG.(names{n}).CS7_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS8_mean=Pwrnan.TPHNEG.(names{n}).CS8_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS9_mean=Pwrnan.TPHNEG.(names{n}).CS9_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHNEG.(names{n}).norm_CS10_mean=Pwrnan.TPHNEG.(names{n}).CS10_mean(:,:)./Pwrnan.TPHNEG.(names{n}).pCS3_10_mean(:,:);

end
end

fooTPHPOS_BNST=LFPsallnan.TPHPOS;
for a=1:length(TPHPOS_BNST_L)
     for n=2%1:length(names)
    fooTPHPOS_BNST=LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{n}).pips.Pips_freqtime;
Pwrnan.TPHPOS.(names{n}).CS1_4_mean(:,a)=fooTPHPOS_BNST.CS1_4_mean;
Pwrnan.TPHPOS.(names{n}).CS7_10_mean(:,a)=fooTPHPOS_BNST.CS7_10_mean;
Pwrnan.TPHPOS.(names{n}).all_CSs_mean(:,a)=fooTPHPOS_BNST.all_CSs_mean;
Pwrnan.TPHPOS.(names{n}).CS1_mean(:,a)=fooTPHPOS_BNST.CS1only;
Pwrnan.TPHPOS.(names{n}).CS2_mean(:,a)=fooTPHPOS_BNST.CS2only;
Pwrnan.TPHPOS.(names{n}).CS3_mean(:,a)=fooTPHPOS_BNST.CS3only;
Pwrnan.TPHPOS.(names{n}).CS4_mean(:,a)=fooTPHPOS_BNST.CS4only;
Pwrnan.TPHPOS.(names{n}).CS5_mean(:,a)=fooTPHPOS_BNST.CS5only;
Pwrnan.TPHPOS.(names{n}).CS6_mean(:,a)=fooTPHPOS_BNST.CS6only;
Pwrnan.TPHPOS.(names{n}).CS7_mean(:,a)=fooTPHPOS_BNST.CS7only;
Pwrnan.TPHPOS.(names{n}).CS8_mean(:,a)=fooTPHPOS_BNST.CS8only;
Pwrnan.TPHPOS.(names{n}).CS9_mean(:,a)=fooTPHPOS_BNST.CS9only;
Pwrnan.TPHPOS.(names{n}).CS10_mean(:,a)=fooTPHPOS_BNST.CS10only;

    fooTPHPOS_BNST=LFPsallnan.TPHPOS.(TPHPOS_BNST_L{a}).TonePwr.(names{n}).pips.Prepips_freqtime; 
Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,a)=fooTPHPOS_BNST.pCS3_10_mean;

Pwrnan.TPHPOS.(names{n}).normCS1_4_mean=Pwrnan.TPHPOS.(names{n}).CS1_4_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).normCS7_10_mean=Pwrnan.TPHPOS.(names{n}).CS7_10_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:); 
Pwrnan.TPHPOS.(names{n}).norm_all_CSs_mean=Pwrnan.TPHPOS.(names{n}).all_CSs_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS1_mean=Pwrnan.TPHPOS.(names{n}).CS1_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS2_mean=Pwrnan.TPHPOS.(names{n}).CS2_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS3_mean=Pwrnan.TPHPOS.(names{n}).CS3_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS4_mean=Pwrnan.TPHPOS.(names{n}).CS4_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS5_mean=Pwrnan.TPHPOS.(names{n}).CS5_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS6_mean=Pwrnan.TPHPOS.(names{n}).CS6_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS7_mean=Pwrnan.TPHPOS.(names{n}).CS7_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS8_mean=Pwrnan.TPHPOS.(names{n}).CS8_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS9_mean=Pwrnan.TPHPOS.(names{n}).CS9_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);
Pwrnan.TPHPOS.(names{n}).norm_CS10_mean=Pwrnan.TPHPOS.(names{n}).CS10_mean(:,:)./Pwrnan.TPHPOS.(names{n}).pCS3_10_mean(:,:);

end
end

%%
%Absolute Power for each region, CS and pCS
what={'all_CSs_mean', 'pCS3_10_mean','CS1_4_mean','CS7_10_mean'};

for n=1:length(names)
    for t=1:length(what)
     Pwrnan_all.(names{n}).TPHPOS.(what{t})=(Pwrnan.TPHPOS.(names{n}).(what{t}));
     Pwrnan_all.(names{n}).TPHNEG.(what{t})=(Pwrnan.TPHNEG.(names{n}).(what{t}));
    end
end

Freq=Freq_pw;
gamma_freq_idx_90_140=find(Freq>=70 & Freq<=150);

figure
n=1 % cea

 subplot(3,4,1)     
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{1})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]); ylim([1 80000]) ;title([char(names(n)), ' RECALL TPH-' ]); 
        
        subplot(3,4,5)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{3})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 80000]); title([char(names(n)), ' RECALL TPH-' ]);
        
        subplot(3,4,9)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{4})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 80000]); title([char(names(n)), ' RECALL TPH-']); 

%TPH+ ACSF
     subplot(3,4,2)     
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{1})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]); ylim([1 80000]) ;title([char(names(n)), ' RECALL TPH+' ]); 
        
        subplot(3,4,6)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{3})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 80000]); title([char(names(n)), ' RECALL TPH+' ]);
        
        subplot(3,4,10)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{4})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 80000]); title([char(names(n)), ' RECALL TPH+']); 
       
        %TPH- ACSF
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        %CeA L
        subplot(3,4,3)
        CS_SEM = std(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH- all CSs ', ' ','TPH- pCS3-10'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,7)
        CS_SEM = std(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
         ylim([0 80000]);legend('TPH- CS1-4 ', ' ','TPH- pCS3-10'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,11)
        CS_SEM = std(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH- CS7-10 ', ' ','TPH- pCS3-10'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
   
        %TPH+ acsf
       fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        %CeA L
        subplot(3,4,4)
        CS_SEM = std(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH+ all CSs ', ' ','TPH+ pCS3-10'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,8)
        CS_SEM = std(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
         ylim([0 80000]);legend('TPH+ CS1-4 ', ' ','TPH+ pCS3-10'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,12)
        CS_SEM = std(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 80000]);legend('TPH+ CS7-10 ', ' ','TPH+ pCS3-10'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])

%%
figure
n=2 % BNST

 subplot(3,4,1)     
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{1})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]); ylim([1 60000]) ;title([char(names(n)), ' RECALL TPH-' ]); 
        
        subplot(3,4,5)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{3})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 60000]); title([char(names(n)), ' RECALL TPH-' ]);
        
        subplot(3,4,9)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{4})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 60000]); title([char(names(n)), ' RECALL TPH-']); 

%TPH+ ACSF
     subplot(3,4,2)     
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{1})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]); ylim([1 60000]) ;title([char(names(n)), ' RECALL TPH+' ]); 
        
        subplot(3,4,6)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{3})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 60000]); title([char(names(n)), ' RECALL TPH+' ]);
        
        subplot(3,4,10)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{4})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([70 150]);ylim([1 60000]); title([char(names(n)), ' RECALL TPH+']); 
       
        %TPH- ACSF
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        %CeA L
        subplot(3,4,3)
        CS_SEM = std(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 50000]);legend('TPH- all CSs ', ' ','TPH- pCS3-10'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,7)
        CS_SEM = std(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
         ylim([0 50000]);legend('TPH- CS1-4 ', ' ','TPH- pCS3-10'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,11)
        CS_SEM = std(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHNEG{4,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 50000]);legend('TPH- CS7-10 ', ' ','TPH- pCS3-10'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
   
        %TPH+ acsf
       fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        %CeA L
        subplot(3,4,4)
        CS_SEM = std(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 50000]);legend('TPH+ all CSs ', ' ','TPH+ pCS3-10'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,8)
        CS_SEM = std(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
         ylim([0 50000]);legend('TPH+ CS1-4 ', ' ','TPH+ pCS3-10'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
        
        subplot(3,4,12)
        CS_SEM = std(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))-1);
        bar(3,(mean(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5469 0')
        hold on; errorbar(3,mean(mean(fooTPHPOS{4,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on; CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);         
        bar(4,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on; errorbar(4,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')     
        ylim([0 50000]);legend('TPH+ CS7-10 ', ' ','TPH+ pCS3-10'); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])

%%
%NORM Power for each region, CS and pCS
clear what
clear Pwrnan_all
what={'norm_all_CSs_mean', 'normCS1_4_mean', 'normCS7_10_mean'};

for n=1:length(names)
    for t=1:length(what)
     Pwrnan_all.(names{n}).TPHPOS.(what{t})=(Pwrnan.TPHPOS.(names{n}).(what{t}));
     Pwrnan_all.(names{n}).TPHNEG.(what{t})=(Pwrnan.TPHNEG.(names{n}).(what{t}));
    end
end

figure(23)   
n=1 ;%cea
        subplot(3,4,1)     
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{1})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{1})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')
       
        subplot(3,4,5)  
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')
        
        subplot(3,4,9)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{3})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{3})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')

subplot(3,4,2)
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')  

        fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1.5]);legend('TPH- all CS ',' ','TPH+ all CS ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
      
        subplot(3,4,6)
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
  
        fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
           ylim([0 1.5]);legend('TPH- CS1-4 ',' ','TPH+ CS1-4 ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
      
        subplot(3,4,10)
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        

        fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1.5]);legend('TPH- CS7-10 ', ' ','TPH+ CS7-10 ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])      

n=2 ;%bnst
        subplot(3,4,3)     
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{1})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{1})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')
       
        subplot(3,4,7)  
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{2})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{2})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')
        
        subplot(3,4,11)
        errorbarplot_joe(Freq_pw,Pwrnan.TPHNEG.(names{n}).(what{3})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_pw,Pwrnan.TPHPOS.(names{n}).(what{3})',[0 0 1],[0 0 1])
        hold on
        xlim([70 150]);ylim([0.5 1.5]); title([char(names(n)), ' RECALL ']); ylabel('NORM Power')

subplot(3,4,4)
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')  

        fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1.5]);legend('TPH- all CS ',' ','TPH+ all CS ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
      
        subplot(3,4,8)
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
  
        fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
           ylim([0 1.5]);legend('TPH- CS1-4 ',' ','TPH+ CS1-4 ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])
      
        subplot(3,4,12)
        fooTPHNEG=struct2cell(Pwrnan_all.(names{n}).TPHNEG);
        CS_SEM = std(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on; errorbar(1,mean(mean(fooTPHNEG{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        
        fooTPHPOS=struct2cell(Pwrnan_all.(names{n}).TPHPOS);
        CS_SEM = std(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on; errorbar(2,mean(mean(fooTPHPOS{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1.5]);legend('TPH- CS7-10 ', ' ','TPH+ CS7-10 ', ' '); 
        title([char(names(n)), ' RECALL Power under curve 90to140Hz '])      

%% Calculating Correlation Between Regions

TPHNEG_BNST_CeA=intersect(TPHNEG_BNST_L,TPHNEG_CeA_L);
TPHPOS_BNST_CeA=intersect(TPHPOS_BNST_L,TPHPOS_CeA_L);

%% merged ctrls Coherence Calculation and restructure
%analyze first for tones on mtchg, then cut for pips
%TPHNEG

for a=1:length(TPHNEG_BNST_CeA)     
          LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Coh.BNSTCeA=struct('twoD_cohe',[]);
          LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Coh.BNSTCeA.twoD_cohe=struct('CS_freqtime',[],'pCS_freqtime',[]);
          LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Coh.BNSTCeA.twoD_cohe.CS_freqtime=struct('CS1',[],'CS2',[],'CS3',[],'CS4',[],'CS5',[],'CS6',[],'CS7',[],'CS8',[],'CS9',[],'CS10',[]);
          LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Coh.BNSTCeA.twoD_cohe.pCS_freqtime=struct('pCS1',[],'pCS2',[],'pCS3',[],'pCS4',[],'pCS5',[],'pCS6',[],'pCS7',[],'pCS8',[],'pCS9',[],'pCS10',[]);
end

BNSTCeA_inputmatrix=struct('CS1',[],'CS2',[],'CS3',[],'CS4',[],'CS5',[],'CS6',[],'CS7',[],'CS8',[],'CS9',[],'CS10',[],'pCS1',[],'pCS2',[],'pCS3',[],'pCS4',[],'pCS5',[],'pCS6',[],'pCS7',[],'pCS8',[],'pCS9',[],'pCS10',[]);

for a=1:length(TPHNEG_BNST_CeA) %CHANGE BELOW FOR REGIONS
    for c=1:numCS              
              BNSTCeA_inputmatrix.(CSname{c}).(TPHNEG_BNST_CeA{a})(1,:)=LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Filt.CeA_L.CS{1, c};
              BNSTCeA_inputmatrix.(CSname{c}).(TPHNEG_BNST_CeA{a})(2,:)=LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Filt.BNST_L.CS{1, c};
              BNSTCeA_inputmatrix.(pCSname{c}).(TPHNEG_BNST_CeA{a})(1,:)=LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Filt.CeA_L.pCS{1, c};
              BNSTCeA_inputmatrix.(pCSname{c}).(TPHNEG_BNST_CeA{a})(2,:)=LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Filt.BNST_L.pCS{1, c};
    end
end

for c=1:numCS
        for a=1:length(TPHNEG_BNST_CeA) 
            [LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Coh.BNSTCeA.CS_freqtime{c},Freq_c,Time_c]=mtchg(BNSTCeA_inputmatrix.(CSname{c}).(TPHNEG_BNST_CeA{a})',2048,2000,500,480,2);  
            [LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Coh.BNSTCeA.pCS_freqtime{c},Freq_c,Time_c]=mtchg(BNSTCeA_inputmatrix.(pCSname{c}).(TPHNEG_BNST_CeA{a})',2048,2000,500,480,2);  
        end
end       

  for a=1:length(TPHNEG_BNST_CeA) 
      for c=1:numCS 
            LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Coh.BNSTCeA.twoD_cohe.CS_freqtime.(CSname{c})=LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Coh.BNSTCeA.CS_freqtime{1,c}(:,:,1,2);
            LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Coh.BNSTCeA.twoD_cohe.pCS_freqtime.(pCSname{c})=LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Coh.BNSTCeA.pCS_freqtime{1,c}(:,:,1,2);
        
      end
  end


%now we isolate the coherence for the bins corresponding to each pip &
%prepip
 for a=1:length(TPHNEG_BNST_CeA)
        for c=1:numCS
            for d=1:num_pips     
                 LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.(CSname{c}).(Pipname{d})=LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Coh.BNSTCeA.twoD_cohe.CS_freqtime.(CSname{c})(:,prepips_bins(d,1):prepips_bins(d,100));
                 LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d})=LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).Coh.BNSTCeA.twoD_cohe.pCS_freqtime.(pCSname{c})(:,prepips_bins(d,1):prepips_bins(d,100));
            end
        end
 end
 
 for a=1:length(TPHNEG_BNST_CeA)
        for c=1:numCS            
            ooooo=LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.(CSname{c});
            LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.(CSname{c}).altogether=cat(2,ooooo.pip1,ooooo.pip2,ooooo.pip3,ooooo.pip4,ooooo.pip5,ooooo.pip6,ooooo.pip7,ooooo.pip8,ooooo.pip9,ooooo.pip10,ooooo.pip11,ooooo.pip12,ooooo.pip13,ooooo.pip14,ooooo.pip15,ooooo.pip16,ooooo.pip17,ooooo.pip18,ooooo.pip19,ooooo.pip20,ooooo.pip21,ooooo.pip22,ooooo.pip23,ooooo.pip24,ooooo.pip25,ooooo.pip26,ooooo.pip27,ooooo.pip28,ooooo.pip29,ooooo.pip30);
        
            ooooo=LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime.(pCSname{c});
            LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime.(pCSname{c}).altogether=cat(2,ooooo.prepip1,ooooo.prepip2,ooooo.prepip3,ooooo.prepip4,ooooo.prepip5,ooooo.prepip6,ooooo.prepip7,ooooo.prepip8,ooooo.prepip9,ooooo.prepip10,ooooo.prepip11,ooooo.prepip12,ooooo.prepip13,ooooo.prepip14,ooooo.prepip15,ooooo.prepip16,ooooo.prepip17,ooooo.prepip18,ooooo.prepip19,ooooo.prepip20,ooooo.prepip21,ooooo.prepip22,ooooo.prepip23,ooooo.prepip24,ooooo.prepip25,ooooo.prepip26,ooooo.prepip27,ooooo.prepip28,ooooo.prepip29,ooooo.prepip30);
       end
 end

for a=1:length(TPHNEG_BNST_CeA)
    ooooo=LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime;
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.all_CSs=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether,ooooo.CS5.altogether,ooooo.CS6.altogether,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.all_CSs_mean=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.all_CSs,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1_4=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1_4_mean=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1_4,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1_2=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1_2_mean=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1_2,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS7_10=cat(2,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS7_10_mean=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS7_10,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS9_10=cat(2,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS9_10_mean=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS9_10,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS2only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS2.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS3only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS3.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS4only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS4.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS5only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS5.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS6only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS6.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS7only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS7.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS8only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS8.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS9only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS9.altogether,2);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS10only=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS10.altogether,2);
end

for a=1:length(TPHNEG_BNST_CeA)
    ooooo=LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime;
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime.all_pCSs=cat(2,ooooo.pCS1.altogether,ooooo.pCS2.altogether,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime.all_pCSs_mean=mean(LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime.all_pCSs,2);
end

%TPHPOS
for a=1:length(TPHPOS_BNST_CeA)     
          LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Coh.BNSTCeA=struct('twoD_cohe',[]);
          LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Coh.BNSTCeA.twoD_cohe=struct('CS_freqtime',[],'pCS_freqtime',[]);
          LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Coh.BNSTCeA.twoD_cohe.CS_freqtime=struct('CS1',[],'CS2',[],'CS3',[],'CS4',[],'CS5',[],'CS6',[],'CS7',[],'CS8',[],'CS9',[],'CS10',[]);
          LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Coh.BNSTCeA.twoD_cohe.pCS_freqtime=struct('pCS1',[],'pCS2',[],'pCS3',[],'pCS4',[],'pCS5',[],'pCS6',[],'pCS7',[],'pCS8',[],'pCS9',[],'pCS10',[]);
end

BNSTCeA_inputmatrix=struct('CS1',[],'CS2',[],'CS3',[],'CS4',[],'CS5',[],'CS6',[],'CS7',[],'CS8',[],'CS9',[],'CS10',[],'pCS1',[],'pCS2',[],'pCS3',[],'pCS4',[],'pCS5',[],'pCS6',[],'pCS7',[],'pCS8',[],'pCS9',[],'pCS10',[]);

for a=1:length(TPHPOS_BNST_CeA) %CHANGE BELOW FOR REGIONS
    for c=1:numCS              
              BNSTCeA_inputmatrix.(CSname{c}).(TPHPOS_BNST_CeA{a})(1,:)=LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Filt.CeA_L.CS{1, c};
              BNSTCeA_inputmatrix.(CSname{c}).(TPHPOS_BNST_CeA{a})(2,:)=LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Filt.BNST_L.CS{1, c};
              BNSTCeA_inputmatrix.(pCSname{c}).(TPHPOS_BNST_CeA{a})(1,:)=LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Filt.CeA_L.pCS{1, c};
              BNSTCeA_inputmatrix.(pCSname{c}).(TPHPOS_BNST_CeA{a})(2,:)=LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Filt.BNST_L.pCS{1, c};
    end
end

for c=1:numCS
        for a=1:length(TPHPOS_BNST_CeA) 
            [LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Coh.BNSTCeA.CS_freqtime{c},Freq_c,Time_c]=mtchg(BNSTCeA_inputmatrix.(CSname{c}).(TPHPOS_BNST_CeA{a})',2048,2000,500,480,2);  
            [LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Coh.BNSTCeA.pCS_freqtime{c},Freq_c,Time_c]=mtchg(BNSTCeA_inputmatrix.(pCSname{c}).(TPHPOS_BNST_CeA{a})',2048,2000,500,480,2);  
        end
end       

  for a=1:length(TPHPOS_BNST_CeA) 
      for c=1:numCS 
            LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Coh.BNSTCeA.twoD_cohe.CS_freqtime.(CSname{c})=LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Coh.BNSTCeA.CS_freqtime{1,c}(:,:,1,2);
            LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Coh.BNSTCeA.twoD_cohe.pCS_freqtime.(pCSname{c})=LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Coh.BNSTCeA.pCS_freqtime{1,c}(:,:,1,2);
        
      end
  end


%now we isolate the coherence for the bins corresponding to each pip &
%prepip
 for a=1:length(TPHPOS_BNST_CeA)
        for c=1:numCS
            for d=1:num_pips     
                 LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.(CSname{c}).(Pipname{d})=LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Coh.BNSTCeA.twoD_cohe.CS_freqtime.(CSname{c})(:,prepips_bins(d,1):prepips_bins(d,100));
                 LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime.(pCSname{c}).(Prepipname{d})=LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).Coh.BNSTCeA.twoD_cohe.pCS_freqtime.(pCSname{c})(:,prepips_bins(d,1):prepips_bins(d,100));
            end
        end
 end
 
 for a=1:length(TPHPOS_BNST_CeA)
        for c=1:numCS            
            ooooo=LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.(CSname{c});
            LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.(CSname{c}).altogether=cat(2,ooooo.pip1,ooooo.pip2,ooooo.pip3,ooooo.pip4,ooooo.pip5,ooooo.pip6,ooooo.pip7,ooooo.pip8,ooooo.pip9,ooooo.pip10,ooooo.pip11,ooooo.pip12,ooooo.pip13,ooooo.pip14,ooooo.pip15,ooooo.pip16,ooooo.pip17,ooooo.pip18,ooooo.pip19,ooooo.pip20,ooooo.pip21,ooooo.pip22,ooooo.pip23,ooooo.pip24,ooooo.pip25,ooooo.pip26,ooooo.pip27,ooooo.pip28,ooooo.pip29,ooooo.pip30);
        
            ooooo=LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime.(pCSname{c});
            LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime.(pCSname{c}).altogether=cat(2,ooooo.prepip1,ooooo.prepip2,ooooo.prepip3,ooooo.prepip4,ooooo.prepip5,ooooo.prepip6,ooooo.prepip7,ooooo.prepip8,ooooo.prepip9,ooooo.prepip10,ooooo.prepip11,ooooo.prepip12,ooooo.prepip13,ooooo.prepip14,ooooo.prepip15,ooooo.prepip16,ooooo.prepip17,ooooo.prepip18,ooooo.prepip19,ooooo.prepip20,ooooo.prepip21,ooooo.prepip22,ooooo.prepip23,ooooo.prepip24,ooooo.prepip25,ooooo.prepip26,ooooo.prepip27,ooooo.prepip28,ooooo.prepip29,ooooo.prepip30);
       end
 end

for a=1:length(TPHPOS_BNST_CeA)
    ooooo=LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime;
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.all_CSs=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether,ooooo.CS5.altogether,ooooo.CS6.altogether,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.all_CSs_mean=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.all_CSs,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1_4=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether,ooooo.CS3.altogether,ooooo.CS4.altogether);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1_4_mean=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1_4,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1_2=cat(2,ooooo.CS1.altogether,ooooo.CS2.altogether);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1_2_mean=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1_2,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS7_10=cat(2,ooooo.CS7.altogether,ooooo.CS8.altogether,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS7_10_mean=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS7_10,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS9_10=cat(2,ooooo.CS9.altogether,ooooo.CS10.altogether);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS9_10_mean=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS9_10,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS1.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS2only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS2.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS3only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS3.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS4only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS4.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS5only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS5.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS6only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS6.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS7only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS7.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS8only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS8.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS9only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS9.altogether,2);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS10only=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime.CS10.altogether,2);
end

for a=1:length(TPHPOS_BNST_CeA)
    ooooo=LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime;
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime.all_pCSs=cat(2,ooooo.pCS1.altogether,ooooo.pCS2.altogether,ooooo.pCS3.altogether,ooooo.pCS4.altogether,ooooo.pCS5.altogether,ooooo.pCS6.altogether,ooooo.pCS7.altogether,ooooo.pCS8.altogether,ooooo.pCS9.altogether,ooooo.pCS10.altogether);
    LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime.all_pCSs_mean=mean(LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime.all_pCSs,2);
end

%% Restructuring Coherence for each group & calculating CS normalized by pcS

for a=1:length(TPHNEG_BNST_CeA)
    fooctrls=LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime;
Coh.TPHNEG.avgs.BNSTCeA.CS1_4_mean(:,a)=fooctrls.CS1_4_mean;
Coh.TPHNEG.avgs.BNSTCeA.CS1_2_mean(:,a)=fooctrls.CS1_2_mean;
Coh.TPHNEG.avgs.BNSTCeA.CS7_10_mean(:,a)=fooctrls.CS7_10_mean; 
Coh.TPHNEG.avgs.BNSTCeA.CS9_10_mean(:,a)=fooctrls.CS9_10_mean; 
Coh.TPHNEG.avgs.BNSTCeA.CS1_mean(:,a)=fooctrls.CS1only;
Coh.TPHNEG.avgs.BNSTCeA.CS2_mean(:,a)=fooctrls.CS2only;
Coh.TPHNEG.avgs.BNSTCeA.CS3_mean(:,a)=fooctrls.CS3only;
Coh.TPHNEG.avgs.BNSTCeA.CS4_mean(:,a)=fooctrls.CS4only;
Coh.TPHNEG.avgs.BNSTCeA.CS5_mean(:,a)=fooctrls.CS5only;
Coh.TPHNEG.avgs.BNSTCeA.CS6_mean(:,a)=fooctrls.CS6only;
Coh.TPHNEG.avgs.BNSTCeA.CS7_mean(:,a)=fooctrls.CS7only;
Coh.TPHNEG.avgs.BNSTCeA.CS8_mean(:,a)=fooctrls.CS8only;
Coh.TPHNEG.avgs.BNSTCeA.CS9_mean(:,a)=fooctrls.CS9only;
Coh.TPHNEG.avgs.BNSTCeA.CS10_mean(:,a)=fooctrls.CS10only;

Coh.TPHNEG.avgs.BNSTCeA.all_CSs_mean(:,a)=fooctrls.all_CSs_mean;

    fooctrls=LFPsall.TPHNEG.(TPHNEG_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime; 
Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,a)=fooctrls.all_pCSs_mean;

Coh.TPHNEG.avgs.BNSTCeA.normCS1_2_mean=Coh.TPHNEG.avgs.BNSTCeA.CS1_2_mean(:,:)./Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHNEG.avgs.BNSTCeA.normCS1_4_mean=Coh.TPHNEG.avgs.BNSTCeA.CS1_4_mean(:,:)./Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHNEG.avgs.BNSTCeA.normCS7_10_mean=Coh.TPHNEG.avgs.BNSTCeA.CS7_10_mean(:,:)./Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,:); 
Coh.TPHNEG.avgs.BNSTCeA.normCS9_10_mean=Coh.TPHNEG.avgs.BNSTCeA.CS9_10_mean(:,:)./Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,:); 
Coh.TPHNEG.avgs.BNSTCeA.norm_all_CSs_mean=Coh.TPHNEG.avgs.BNSTCeA.all_CSs_mean(:,:)./Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHNEG.avgs.BNSTCeA.normCS1_mean=Coh.TPHNEG.avgs.BNSTCeA.CS1_mean(:,:)./Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHNEG.avgs.BNSTCeA.normCS2_mean=Coh.TPHNEG.avgs.BNSTCeA.CS2_mean(:,:)./Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHNEG.avgs.BNSTCeA.normCS3_mean=Coh.TPHNEG.avgs.BNSTCeA.CS3_mean(:,:)./Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHNEG.avgs.BNSTCeA.normCS4_mean=Coh.TPHNEG.avgs.BNSTCeA.CS4_mean(:,:)./Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHNEG.avgs.BNSTCeA.normCS5_mean=Coh.TPHNEG.avgs.BNSTCeA.CS5_mean(:,:)./Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHNEG.avgs.BNSTCeA.normCS6_mean=Coh.TPHNEG.avgs.BNSTCeA.CS6_mean(:,:)./Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHNEG.avgs.BNSTCeA.normCS7_mean=Coh.TPHNEG.avgs.BNSTCeA.CS7_mean(:,:)./Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHNEG.avgs.BNSTCeA.normCS8_mean=Coh.TPHNEG.avgs.BNSTCeA.CS8_mean(:,:)./Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHNEG.avgs.BNSTCeA.normCS9_mean=Coh.TPHNEG.avgs.BNSTCeA.CS9_mean(:,:)./Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHNEG.avgs.BNSTCeA.normCS10_mean=Coh.TPHNEG.avgs.BNSTCeA.CS10_mean(:,:)./Coh.TPHNEG.avgs.BNSTCeA.all_pCSs_mean(:,:);

end

%TPHPOS
for a=1:length(TPHPOS_BNST_CeA)
    fooctrls=LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Pips_freqtime;
Coh.TPHPOS.avgs.BNSTCeA.CS1_4_mean(:,a)=fooctrls.CS1_4_mean;
Coh.TPHPOS.avgs.BNSTCeA.CS1_2_mean(:,a)=fooctrls.CS1_2_mean;
Coh.TPHPOS.avgs.BNSTCeA.CS7_10_mean(:,a)=fooctrls.CS7_10_mean; 
Coh.TPHPOS.avgs.BNSTCeA.CS9_10_mean(:,a)=fooctrls.CS9_10_mean; 
Coh.TPHPOS.avgs.BNSTCeA.CS1_mean(:,a)=fooctrls.CS1only;
Coh.TPHPOS.avgs.BNSTCeA.CS2_mean(:,a)=fooctrls.CS2only;
Coh.TPHPOS.avgs.BNSTCeA.CS3_mean(:,a)=fooctrls.CS3only;
Coh.TPHPOS.avgs.BNSTCeA.CS4_mean(:,a)=fooctrls.CS4only;
Coh.TPHPOS.avgs.BNSTCeA.CS5_mean(:,a)=fooctrls.CS5only;
Coh.TPHPOS.avgs.BNSTCeA.CS6_mean(:,a)=fooctrls.CS6only;
Coh.TPHPOS.avgs.BNSTCeA.CS7_mean(:,a)=fooctrls.CS7only;
Coh.TPHPOS.avgs.BNSTCeA.CS8_mean(:,a)=fooctrls.CS8only;
Coh.TPHPOS.avgs.BNSTCeA.CS9_mean(:,a)=fooctrls.CS9only;
Coh.TPHPOS.avgs.BNSTCeA.CS10_mean(:,a)=fooctrls.CS10only;

Coh.TPHPOS.avgs.BNSTCeA.all_CSs_mean(:,a)=fooctrls.all_CSs_mean;

    fooctrls=LFPsall.TPHPOS.(TPHPOS_BNST_CeA{a}).ToneCoh.BNSTCeA.pips.Prepips_freqtime; 
Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,a)=fooctrls.all_pCSs_mean;

Coh.TPHPOS.avgs.BNSTCeA.normCS1_2_mean=Coh.TPHPOS.avgs.BNSTCeA.CS1_2_mean(:,:)./Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHPOS.avgs.BNSTCeA.normCS1_4_mean=Coh.TPHPOS.avgs.BNSTCeA.CS1_4_mean(:,:)./Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHPOS.avgs.BNSTCeA.normCS7_10_mean=Coh.TPHPOS.avgs.BNSTCeA.CS7_10_mean(:,:)./Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,:); 
Coh.TPHPOS.avgs.BNSTCeA.normCS9_10_mean=Coh.TPHPOS.avgs.BNSTCeA.CS9_10_mean(:,:)./Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,:); 
Coh.TPHPOS.avgs.BNSTCeA.norm_all_CSs_mean=Coh.TPHPOS.avgs.BNSTCeA.all_CSs_mean(:,:)./Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHPOS.avgs.BNSTCeA.normCS1_mean=Coh.TPHPOS.avgs.BNSTCeA.CS1_mean(:,:)./Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHPOS.avgs.BNSTCeA.normCS2_mean=Coh.TPHPOS.avgs.BNSTCeA.CS2_mean(:,:)./Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHPOS.avgs.BNSTCeA.normCS3_mean=Coh.TPHPOS.avgs.BNSTCeA.CS3_mean(:,:)./Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHPOS.avgs.BNSTCeA.normCS4_mean=Coh.TPHPOS.avgs.BNSTCeA.CS4_mean(:,:)./Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHPOS.avgs.BNSTCeA.normCS5_mean=Coh.TPHPOS.avgs.BNSTCeA.CS5_mean(:,:)./Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHPOS.avgs.BNSTCeA.normCS6_mean=Coh.TPHPOS.avgs.BNSTCeA.CS6_mean(:,:)./Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHPOS.avgs.BNSTCeA.normCS7_mean=Coh.TPHPOS.avgs.BNSTCeA.CS7_mean(:,:)./Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHPOS.avgs.BNSTCeA.normCS8_mean=Coh.TPHPOS.avgs.BNSTCeA.CS8_mean(:,:)./Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHPOS.avgs.BNSTCeA.normCS9_mean=Coh.TPHPOS.avgs.BNSTCeA.CS9_mean(:,:)./Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,:);
Coh.TPHPOS.avgs.BNSTCeA.normCS10_mean=Coh.TPHPOS.avgs.BNSTCeA.CS10_mean(:,:)./Coh.TPHPOS.avgs.BNSTCeA.all_pCSs_mean(:,:);

end

%%
%non-norm Coh for each region, CS and pCS

clear Coh_all
clear what

what={'all_CSs_mean', 'all_pCSs_mean','CS1_4_mean','CS7_10_mean'};

names={'BNSTCeA'}

for t=1:length(what)
     Coh_all.BNSTCeA.TPHNEG.(what{t})=(Coh.TPHNEG.avgs.BNSTCeA.(what{t}));
     Coh_all.BNSTCeA.TPHPOS.(what{t})=(Coh.TPHPOS.avgs.BNSTCeA.(what{t}));
end

%change frequency assessed by inserting frequency as indicated below
Freq=Freq_c;

gamma_freq_idx_90_140=find(Freq>=90 & Freq<=140);

for n=1:length(names)
figure(2)     
     subplot(2,3,1)
        hold on
        errorbarplot_joe(Freq_c,Coh.TPHNEG.avgs.BNSTCeA.(what{1})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_c,Coh.TPHNEG.avgs.BNSTCeA.(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([95 150]);ylim([0.4 1]);
        title([char(names(n)), ' TPHNEG ' char(what(1)), ' ' char(what(2))])
        legend('CTRLS all CSs ',' ','CTRLS all pCSs ',' ')
        ylabel('BNST-CeA Coherence ')
        
        subplot(2,3,2)
        errorbarplot_joe(Freq_c,Coh.TPHNEG.avgs.BNSTCeA.(what{3})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_c,Coh.TPHNEG.avgs.BNSTCeA.(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([95 150]);ylim([0.4 1]);
        title([char(names(n)), ' TPHNEG ' char(what(3)), ' ' char(what(2))])
        legend('CS1-4 ',' ','all pCSs ',' ')
        ylabel('BNST-CeA Coherence ')

       subplot(2,3,3)
         errorbarplot_joe(Freq_c,Coh.TPHNEG.avgs.BNSTCeA.(what{4})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_c,Coh.TPHNEG.avgs.BNSTCeA.(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])        
        xlim([95 150]);ylim([0.4 1]);
        title([char(names(n)), ' TPHNEG ' char(what(4)), ' ' char(what(2))])
        legend('CS7-10 ',' ','all pCSs ',' ')
        ylabel('BNST-CeA Coherence ')
       
        %AUC
        fooctrls=struct2cell(Coh_all.BNSTCeA.TPHNEG);
       subplot(2,3,4)
        CS_SEM = std(mean(fooctrls{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooctrls{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5 0')
        hold on 
        errorbar(1,mean(mean(fooctrls{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on
         CS_SEM = std(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{2,1}(gamma_freq_idx_90_140,:)))-1);
         bar(2,(mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
         hold on 
         errorbar(2,mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1]); 
        title([char(names(n)), ' TPHNEG ' char(what(1))])
        legend('all CSs ', ' ','pCSs')

        subplot(2,3,5)
        CS_SEM = std(mean(fooctrls{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooctrls{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5 0')
        hold on 
        errorbar(1,mean(mean(fooctrls{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on
         CS_SEM = std(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{2,1}(gamma_freq_idx_90_140,:)))-1);
         bar(2,(mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
         hold on 
         errorbar(2,mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
          ylim([0 1]); 
        title([char(names(n)), ' TPHNEG ' char(what(3))])
        legend('CS1-4 ', ' ','all pCSs')
         
        subplot(2,3,6)
        CS_SEM = std(mean(fooctrls{4,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{4,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooctrls{4,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5 0')
        hold on 
        errorbar(1,mean(mean(fooctrls{4,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on
        CS_SEM = std(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on 
        errorbar(2,mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1]); 
        title([char(names(n)), ' TPHNEG ' char(what(4))])
        legend('C7-10 ', ' ','all pCSs')
        
end

%tphPOS

for n=1:length(names)
figure(3)     
     subplot(2,3,1)
        hold on
        errorbarplot_joe(Freq_c,Coh.TPHPOS.avgs.BNSTCeA.(what{1})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_c,Coh.TPHPOS.avgs.BNSTCeA.(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([95 150]);ylim([0.4 1]);
        title([char(names(n)), ' TPHPOS ' char(what(1)), ' ' char(what(2))])
        legend('CTRLS all CSs ',' ','CTRLS all pCSs ',' ')
        ylabel('BNST-CeA Coherence ')
        
        subplot(2,3,2)
        errorbarplot_joe(Freq_c,Coh.TPHPOS.avgs.BNSTCeA.(what{3})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_c,Coh.TPHPOS.avgs.BNSTCeA.(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])
        xlim([95 150]);ylim([0.4 1]);
        title([char(names(n)), ' TPHPOS ' char(what(3)), ' ' char(what(2))])
        legend('CS1-4 ',' ','all pCSs ',' ')
        ylabel('BNST-CeA Coherence ')

       subplot(2,3,3)
         errorbarplot_joe(Freq_c,Coh.TPHPOS.avgs.BNSTCeA.(what{4})',[1.0000    0.5469         0],[1.0000    0.5469         0])
        hold on
        errorbarplot_joe(Freq_c,Coh.TPHPOS.avgs.BNSTCeA.(what{2})',[0.5000         0    0.5000],[0.5000         0    0.5000])        
        xlim([95 150]);ylim([0.4 1]);
        title([char(names(n)), ' TPHPOS ' char(what(4)), ' ' char(what(2))])
        legend('CS7-10 ',' ','all pCSs ',' ')
        ylabel('BNST-CeA Coherence ')
       
        %AUC
        fooctrls=struct2cell(Coh_all.BNSTCeA.TPHPOS);
       subplot(2,3,4)
        CS_SEM = std(mean(fooctrls{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooctrls{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5 0')
        hold on 
        errorbar(1,mean(mean(fooctrls{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on
         CS_SEM = std(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{2,1}(gamma_freq_idx_90_140,:)))-1);
         bar(2,(mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
         hold on 
         errorbar(2,mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1]); 
        title([char(names(n)), ' TPHPOS  ' char(what(1))])
        legend('all CSs ', ' ','pCSs')

        subplot(2,3,5)
        CS_SEM = std(mean(fooctrls{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooctrls{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5 0')
        hold on 
        errorbar(1,mean(mean(fooctrls{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on
         CS_SEM = std(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{2,1}(gamma_freq_idx_90_140,:)))-1);
         bar(2,(mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
         hold on 
         errorbar(2,mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
          ylim([0 1]); 
        title([char(names(n)), ' TPHPOS  ' char(what(3))])
        legend('CS1-4 ', ' ','all pCSs')
         
        subplot(2,3,6)
        CS_SEM = std(mean(fooctrls{4,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{4,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooctrls{4,1}(gamma_freq_idx_90_140,:)))),'FaceColor','1.0000 0.5 0')
        hold on 
        errorbar(1,mean(mean(fooctrls{4,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on
        CS_SEM = std(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.5000 0 0.5000')
        hold on 
        errorbar(2,mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        ylim([0 1]); 
        title([char(names(n)), ' TPHPOS ' char(what(4))])
        legend('C7-10 ', ' ','all pCSs')
        
end
%%
%normalized coherence

clear what
clear Coh_all

what={'norm_all_CSs_mean', 'normCS1_4_mean', 'normCS7_10_mean'};

for t=1:length(what)
     Coh_all.BNSTCeA.TPHNEG.(what{t})=(Coh.TPHNEG.avgs.BNSTCeA.(what{t}));
     Coh_all.BNSTCeA.TPHPOS.(what{t})=(Coh.TPHPOS.avgs.BNSTCeA.(what{t}));

end


figure(4)     
        subplot(2,3,1)
        errorbarplot_joe(Freq_c,Coh.TPHNEG.avgs.BNSTCeA.(what{1})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_c,Coh.TPHPOS.avgs.BNSTCeA.(what{1})',[0 0 1],[0 0 1])
        xlim([70 150]);ylim([0.5 1.5]);
        title([char(names(n)), ' NORM COHERENCE '])
        legend('TPH- all CS','','TPH+ all CS','')
        ylabel('BNST-CeA NORM Coherence ')

        subplot(2,3,2)
        errorbarplot_joe(Freq_c,Coh.TPHNEG.avgs.BNSTCeA.(what{2})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_c,Coh.TPHPOS.avgs.BNSTCeA.(what{2})',[0 0 1],[0 0 1])
        xlim([70 150]);ylim([0.5 1.5]);
        title([char(names(n)), ' NORM COHERENCE '])
        legend(['TPH- CS1-4'],'','TPH+ CS1-4','')
        ylabel('BNST-CeA NORM Coherence ')
    
        subplot(2,3,3)
        errorbarplot_joe(Freq_c,Coh.TPHNEG.avgs.BNSTCeA.(what{3})',[0.8 0.8 0.8],[0.5 0.5 0.5])
        hold on
        errorbarplot_joe(Freq_c,Coh.TPHPOS.avgs.BNSTCeA.(what{3})',[0 0 1],[0 0 1])
        xlim([70 150]);ylim([0.5 1.5]);
        title([char(names(n)), ' NORM COHERENCE '])
        legend(['TPH- CS7-10'],'','TPH+ CS7-10','')
        ylabel('BNST-CeA NORM Coherence ')

        subplot(2,3,4)
        fooctrls=struct2cell(Coh_all.BNSTCeA.TPHNEG);
        CS_SEM = std(mean(fooctrls{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooctrls{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on 
        errorbar(1,mean(mean(fooctrls{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on
        fooctrls=struct2cell(Coh_all.BNSTCeA.TPHPOS);
        CS_SEM = std(mean(fooctrls{1,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{1,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooctrls{1,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on 
        errorbar(2,mean(mean(fooctrls{1,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on
        ylim([0 1.2]); 
        title([char(names(n)), ' NORM COHERENCE  ' char(what(1))])
        legend('TPH- all CS','','TPHPOS all CS')

        subplot(2,3,5)
        fooctrls=struct2cell(Coh_all.BNSTCeA.TPHNEG);
        CS_SEM = std(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on 
        errorbar(1,mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on
        fooctrls=struct2cell(Coh_all.BNSTCeA.TPHPOS);
        CS_SEM = std(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{2,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on 
        errorbar(2,mean(mean(fooctrls{2,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on
        ylim([0 1.2]); 
        title([char(names(n)), ' NORM COHERENCE  ' char(what(2))])
        legend('TPH- all CS','','TPHPOS all CS')

        subplot(2,3,6)
        fooctrls=struct2cell(Coh_all.BNSTCeA.TPHNEG);
        CS_SEM = std(mean(fooctrls{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(1,(mean(mean(fooctrls{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0.8 0.8 0.8')
        hold on 
        errorbar(1,mean(mean(fooctrls{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on
        fooctrls=struct2cell(Coh_all.BNSTCeA.TPHPOS);
        CS_SEM = std(mean(fooctrls{3,1}(gamma_freq_idx_90_140,:)))/sqrt(length(sum(fooctrls{3,1}(gamma_freq_idx_90_140,:)))-1);
        bar(2,(mean(mean(fooctrls{3,1}(gamma_freq_idx_90_140,:)))),'FaceColor','0 0 1')
        hold on 
        errorbar(2,mean(mean(fooctrls{3,1}(gamma_freq_idx_90_140,:))),mean(CS_SEM), 'k.')
        hold on
        ylim([0 1.2]); 
        title([char(names(n)), ' NORM COHERENCE  ' char(what(3))])
        legend('TPH- all CS','','TPHPOS all CS')
