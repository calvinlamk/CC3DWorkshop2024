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

#YG Signaling Parameters
ALPHAYG=1 
BETAYG=1750 
EPSILONYG=1000 
KAPPAYG=25000 
THRESHOLDUPYG=5000 
THRESHOLDDOYG=5000

#BR Signaling Parameters
ALPHABR=1 
BETABR=1000 
EPSILONBR=500
KAPPABR=25000 
THRESHOLDUPBR=5000 
THRESHOLDDOBR=5000

#Sampling and Computational Threads Parameters
RESOL=100 #Data sampling rate

class NewSimulationSteppable(SteppableBasePy):

    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)
    def start(self):

#Calling Adhesion Values From the XML File Code   
        global YtoY,YtoG,GtoY,YtoB,BtoY,YtoR,RtoY,GtoG,GtoB,BtoG,GtoR,RtoG,BtoB,BtoR,RtoB,RtoR #the adhesion matrix, call these values and store them for motility code
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
        
#Initialization of Cells Parameters
        for cell in self.cell_list:
            cell.dict["RDM"]=RNG.gauss(RADAVG,RADDEV) #assign the cells a random target radius
            cell.lambdaSurface=2.5                    #temporary value, will be changed in later code
            cell.targetSurface=4*math.pi*cell.dict["RDM"]**2  #spherical surface area
            cell.lambdaVolume=2.5                     #temporary value, will be changed in later code
            cell.targetVolume=(4/3)*math.pi*cell.dict["RDM"]**3 #spherical volume
            cell.dict["PTS"]=[0]                      #initial points for cell, feel free to randomize if desired

    def step(self, mcs):
        
#Defining Per Step Variables Required for Quantifications Code

        for cell in self.cell_list: #iterate over cell list
            CSAY=0 #each cell detect how much surface area it shares with Y cells
            CSAG=0 #each cell detect how much surface area it shares with G cells
            CSAB=0 #each cell detect how much surface area it shares with B cells
            CSAR=0 #each cell detect how much surface area it shares with R cells
            CSAM=0 #each cell detect how much surface area it shares with medium
            
            PTSY=0 #each cell gains points from neighbor type Y
            PTSG=0 #each cell gains points from neighbor type G
            PTSB=0 #each cell gains points from neighbor type B
            PTSR=0 #each cell gains points from neighbor type R
            DTRES=0 #change in reporter due to signal S

            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell): #iterate for each cell its neighbors
                if neighbor is None: #none type refers to medium 
                    continue
                if neighbor.type==1: #gray cells
                    CSAY+=common_surface_area #total common surface area with gray cells
                    PTSY+=0                   #points gray cells send to receiver
                if neighbor.type==2: #green cells
                    CSAG+=common_surface_area #total common surface area with green cells
                    PTSG+=common_surface_area*neighbor.dict["PTS"][0]/(neighbor.surface) #PHI and L as in the text
                if neighbor.type==3: #blue cells
                    CSAB+=common_surface_area #total common surface area with blue cells                     
                    PTSB+=common_surface_area*CONEXPSCF/neighbor.surface
                if neighbor.type==4: #red cells
                    CSAR+=common_surface_area #total common surface area with red cells 
                    PTSR+=common_surface_area*CONEXPSCF/neighbor.surface
            CSAM=cell.surface-(CSAY+CSAG+CSAB+CSAR) #alternative method to calculate common surface area with medium   
            
#Changing Reporter as a Result of Signaling Changes 
            if (cell.type==1 or cell.type==2): #which cells receive what type of signal
                DTRES=(1/(ALPHAYG+math.exp(-((PTSR+PTSB)-BETAYG)/EPSILONYG)))-(1/KAPPAYG)*cell.dict["PTS"][0] #del reporter
                cell.dict["PTS"][0]+=DTRES #update total reporter
#            if (cell.type==3 or cell.type==4): #which cells receive what type of signal
#                DTRES=(1/(ALPHABR+math.exp(-((PTSY+PTSG)-BETABR)/EPSILONBR)))-(1/KAPPABR)*cell.dict["PTS"][0] #del reporter
#                cell.dict["PTS"][0]+=DTRES #update total reporter
                
#Changing State as a Result of Reporter Changes                
            if cell.type==1: #change Y cell state
                if cell.dict["PTS"][0]>=THRESHOLDUPYG:
                    cell.type=2
            if cell.type==2: #change G cell state
                if cell.dict["PTS"][0]<THRESHOLDDOYG:
                    cell.type=1
            if cell.type==3: #change B cell state
                if cell.dict["PTS"][0]>=THRESHOLDUPBR:
                    cell.type=4
            if cell.type==4: #change R cell state
                if cell.dict["PTS"][0]<THRESHOLDDOBR:
                    cell.type=3
                    
#Defining and Calculating Cell Physical Properties Code/Parameters
            if cell.type==1: #gray cells             
                cell.lambdaSurface=2.5            #change depending on cell adhesitivity
                cell.lambdaVolume=2.5             #change depending on cell adhesitivity
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+YtoY*CSAY+YtoG*CSAG+YtoB*CSAB+YtoR*CSAR)/cell.surface #corrected cell motility, tune based on adhesive neighbors, vetted

            if cell.type==2: #green cells
                cell.lambdaSurface=2.5            #change depending on cell adhesitivity
                cell.lambdaVolume=2.5             #change depending on cell adhesitivity    
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+GtoY*CSAY+GtoG*CSAG+GtoB*CSAB+GtoR*CSAR)/cell.surface #corrected cell motility, tune based on adhesive neighbors, vetted
            
            if cell.type==3: #blue cells                              
                cell.lambdaSurface=2.5           #change depending on cell adhesitivity
                cell.lambdaVolume=2.5            #change depending on cell adhesitivity      
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+BtoY*CSAY+BtoG*CSAG+BtoB*CSAB+BtoR*CSAR)/cell.surface # corrected cell motility, vetted
               
            if cell.type==4: #red cells                    
                cell.lambdaSurface=2.5            #change depending on cell adhesitivity
                cell.lambdaVolume=2.5             #change depending on cell adhesitivity      
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+RtoY*CSAY+RtoG*CSAG+RtoB*CSAB+RtoR*CSAR)/cell.surface #corrected cell motility, vetted
              
    def finish(self):
        pass