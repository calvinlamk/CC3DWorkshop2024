from cc3d.core.PySteppables import *
import math                      
import numpy                   
import sys
import random
RNG=random.SystemRandom()

#Global Parameters
#Motility Variables
CtoM=52  #cell to adhesion value
BASAL=100 #baseline motility
SCF=0.5 #scales change in motility

#EndTime Parameters
ENDMCS=50000

#Cell Size and Division Parameters
RADAVG=3 #average radius of the gaussian distribution to choose random radius
RADDEV=.1 #standard deviation of target radius, too low and division couples, too high and you'll lose cells

#Signaling Parameters

#Constitutive Ligand Parameters
CONEXPSCF=10000 #Steady state expression of ligand expressed on a sender cell. This ligand is unaffected by signaling.

#YGR Signaling Parameters
ALPHAYGR=1 
BETAYGR=1750 
EPSILONYGR=1000 
KAPPAYGR=25000 
THRESHOLD=5000 


#BR Signaling Parameters
#ALPHABR=1 
#BETABR=1000 
#EPSILONBR=500
#KAPPABR=25000 
#THRESHOLDUPBR=5000 
#THRESHOLDDOBR=5000

#Sampling and Computational Threads Parameters
RESOL=100 #Data sampling rate

class NewSimulationSteppable(SteppableBasePy):

    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)
    def start(self):

#Calling Adhesion Values From the XML File Code   
        global YtoY,YtoG,GtoY,YtoB,BtoY,YtoR,RtoY,GtoG,GtoB,BtoG,GtoR,RtoG,BtoB,BtoR,RtoB,RtoR,OtoY,YtoO,OtoG,GtoO,OtoB,BtoO,OtoR,RtoO,OtoO 
        YtoYC1=self.get_xml_element('YtoY')
        YtoGC1=GtoYC1=self.get_xml_element('YtoG')
        YtoBC1=BtoYC1=self.get_xml_element('YtoB')
        YtoRC1=RtoYC1=self.get_xml_element('YtoR')
        GtoGC1=self.get_xml_element('GtoG')
        GtoBC1=BtoGC1=self.get_xml_element('GtoB')
        GtoRC1=RtoGC1=self.get_xml_element('GtoR')
        BtoBC1=self.get_xml_element('BtoB')
        BtoRC1=RtoBC1=self.get_xml_element('BtoR')
        RtoRC1=self.get_xml_element('RtoR')
        OtoYC1=YtoOC1=self.get_xml_element('OtoY')
        OtoGC1=GtoOC1=self.get_xml_element('OtoG')
        OtoBC1=BtoOC1=self.get_xml_element('OtoB')
        OtoRC1=RtoOC1=self.get_xml_element('OtoR')
        OtoOC1=self.get_xml_element('OtoO')
        
        YtoY=float(YtoYC1.cdata)
        YtoG=GtoY=float(YtoGC1.cdata)
        YtoB=BtoY=float(YtoBC1.cdata)
        YtoR=RtoY=float(YtoRC1.cdata)
        GtoG=float(GtoGC1.cdata)
        GtoB=BtoG=float(GtoBC1.cdata)
        GtoR=RtoG=float(GtoRC1.cdata)
        BtoB=float(BtoBC1.cdata)
        BtoR=RtoB=float(BtoRC1.cdata)
        RtoR=float(RtoRC1.cdata)
        OtoY=YtoO=float(OtoYC1.cdata)
        OtoG=GtoO=float(OtoGC1.cdata)
        OtoB=BtoO=float(OtoBC1.cdata)
        OtoR=RtoO=float(OtoRC1.cdata)
        OtoO=float(OtoOC1.cdata)  
#Initialization of Cells Parameters
        for cell in self.cell_list:
            cell.dict["RDM"]=RNG.gauss(RADAVG,RADDEV) #assign the cells a random target radius
            cell.lambdaSurface=2.5                    #temporary value, will be changed in later code
            cell.targetSurface=4*math.pi*cell.dict["RDM"]**2  #spherical surface area
            cell.lambdaVolume=2.5                     #temporary value, will be changed in later code
            cell.targetVolume=(4/3)*math.pi*cell.dict["RDM"]**3 #spherical volume
            cell.dict["PTS"]=[0,0]                      #0 for amplifier (red), 1 for adhesion protein

    def step(self, mcs):
        
#Defining Per Step Variables Required for Quantifications Code

        for cell in self.cell_list: #iterate over cell list
            CSAY=0 #each cell detect how much surface area it shares with Y cells
            CSAG=0 #each cell detect how much surface area it shares with G cells
            CSAB=0 #each cell detect how much surface area it shares with B cells
            CSAR=0 #each cell detect how much surface area it shares with R cells
            CSAM=0 #each cell detect how much surface area it shares with medium
            CSAO=0 #each cell detect how much surface area it shares with O cells
            
            PTSY=0 #each cell gains points from neighbor type Y
            PTSG=0 #each cell gains points from neighbor type G
            PTSB=0 #each cell gains points from neighbor type B
            PTSR=0 #each cell gains points from neighbor type R
            PTSO=0 #each cell gains points from neighbor type O
            DTREPO=0 #delta reporter One
            DTREPT=0 #delta reporter Two

            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell): #iterate for each cell its neighbors
                if neighbor is None: #none type refers to medium 
                    continue
                if neighbor.type==1: #gray cells
                    CSAY+=common_surface_area #total common surface area with gray cells
                    PTSY+=0                   #points gray cells send to receiver
                if neighbor.type==2: #green cells
                    CSAG+=common_surface_area #total common surface area with green cells
                    PTSG+=0
                if neighbor.type==3: #blue cells
                    CSAB+=common_surface_area #total common surface area with blue cells                     
                    PTSB+=common_surface_area*CONEXPSCF/neighbor.surface
                if neighbor.type==4: #red cells
                    CSAR+=common_surface_area #total common surface area with red cells 
                    PTSR+=0
                if neighbor.type==5: #orange cells
                    CSAO+=common_surface_area #total common surface area with orange cells 
                    PTSO+=0
            CSAM=cell.surface-(CSAY+CSAG+CSAB+CSAR+CSAO) #alternative method to calculate common surface area with medium   
            
#Changing Reporter as a Result of Signaling Changes 
            if (cell.type==1 or cell.type==2 or cell.type==4 or cell.type==5):
                DTREPO=(1/(ALPHAYGR+math.exp(-((PTSB)-BETAYGR)/EPSILONYGR)))-(1/KAPPAYGR)*cell.dict["PTS"][0] #red amplifier pts
                cell.dict["PTS"][0]+=DTREPO
                DTREPT=(1/(ALPHAYGR+math.exp(-((cell.dict["PTS"][0])-BETAYGR)/EPSILONYGR)))-(1/KAPPAYGR)*cell.dict["PTS"][1] #green adhesion pts
                cell.dict["PTS"][1]+=DTREPT
                
#Changing State as a Result of Reporter Changes                
            if (cell.type==1 or cell.type==2 or cell.type==4 or cell.type==5): #FACS PLOTS AND ADHESION LINKAGE TO STATE
                if cell.dict["PTS"][0]>=THRESHOLD and cell.dict["PTS"][1]<THRESHOLD: #Amplified but not adhesive,r4   
                    cell.type=4
                if cell.dict["PTS"][1]>=THRESHOLD and cell.dict["PTS"][0]<THRESHOLD: #adhesive but not amplifed, g2   
                    cell.type=2                 
                if cell.dict["PTS"][0]>=THRESHOLD and cell.dict["PTS"][1]>=THRESHOLD: #adhesive and amplified, o5  
                    cell.type=5
                if cell.dict["PTS"][0]<THRESHOLD and cell.dict["PTS"][1]<THRESHOLD:   #neither adhesive nor amplifed, g1
                    cell.type=1
                    
#Defining and Calculating Cell Physical Properties Code/Parameters
            if cell.type==1: #Y cells             
                cell.lambdaSurface=2.5            #change depending on cell adhesitivity
                cell.lambdaVolume=2.5             #change depending on cell adhesitivity
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+YtoY*CSAY+YtoG*CSAG+YtoB*CSAB+YtoR*CSAR+YtoO*CSAO)/cell.surface                

            if cell.type==2: #G cells
                cell.lambdaSurface=2.5            #change depending on cell adhesitivity
                cell.lambdaVolume=2.5             #change depending on cell adhesitivity    
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+GtoY*CSAY+GtoG*CSAG+GtoB*CSAB+GtoR*CSAR+GtoO*CSAO)/cell.surface
            
            if cell.type==3: #B cells                              
                cell.lambdaSurface=2.5           #change depending on cell adhesitivity
                cell.lambdaVolume=2.5            #change depending on cell adhesitivity      
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+BtoY*CSAY+BtoG*CSAG+BtoB*CSAB+BtoR*CSAR+BtoO*CSAO)/cell.surface
               
            if cell.type==4: #R cells                    
                cell.lambdaSurface=2.5            #change depending on cell adhesitivity
                cell.lambdaVolume=2.5             #change depending on cell adhesitivity      
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+RtoY*CSAY+RtoG*CSAG+RtoB*CSAB+RtoR*CSAR+RtoO*CSAO)/cell.surface
                
            if cell.type==5: #O cells                    
                cell.lambdaSurface=2.5            #change depending on cell adhesitivity
                cell.lambdaVolume=2.5             #change depending on cell adhesitivity      
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+OtoY*CSAY+OtoG*CSAG+OtoB*CSAB+OtoR*CSAR+OtoO*CSAO)/cell.surface            
              
    def finish(self):
        pass