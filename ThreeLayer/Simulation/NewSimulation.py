from cc3d import CompuCellSetup
        
from .NewSimulationSteppables import NewSimulationSteppable
CompuCellSetup.register_steppable(steppable=NewSimulationSteppable(frequency=1))

CompuCellSetup.run()
