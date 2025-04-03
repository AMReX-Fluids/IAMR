import os
import subprocess

def test_lidDrivenCanvity(working_dir, print_output):
    subprocess.run(
    "make -j8 && ./amr2d.gnu.ex inputs.2d.lid_driven_cavity", 
    shell=True, 
    check=True, 
    cwd=working_dir,
    stdout=None if print_output else subprocess.DEVNULL,
    stderr=None if print_output else subprocess.DEVNULL
    )
    print("Test_lidDrivenCanvity succceed")

def test_LidDrivenCavitySphere(working_dir, print_output):
    subprocess.run(
    "make -j8 && ./amr3d.gnu.MPI.ex inputs.3d.lid_driven_cavity_particle", 
    shell=True, 
    check=True, 
    cwd=working_dir,
    stdout=None if print_output else subprocess.DEVNULL,
    stderr=None if print_output else subprocess.DEVNULL
    )
    print("Test LidDrivenCavitySphere succceed")
    
def test_RSV(working_dir, print_output):
    subprocess.run(
    "make -j8 && ./amr2d.gnu.MPI.ex inputs.2d.rsv", 
    shell=True, 
    check=True, 
    cwd=working_dir,
    stdout=None if print_output else subprocess.DEVNULL,
    stderr=None if print_output else subprocess.DEVNULL
    )
    print("Test RSV succceed")
    
def test_DraftingKissingTumbling(working_dir, print_output):
    subprocess.run(
    "make -j8 USE_CUDA=FALSE USE_MPI=TRUE DEBUG=FALSE;"
    "mpiexec -np 2 ./amr3d.gnu.MPI.ex inputs.3d.DKT max_step=1 amr.n_cell=16 8 8",     
    shell=True, 
    check=True, 
    cwd=working_dir,
    stdout=None if print_output else subprocess.DEVNULL,
    stderr=None if print_output else subprocess.DEVNULL
    )
    print("test_DraftingKissingTumbling succceed")

def test_FallingSphere(working_dir, print_output):
    subprocess.run(
    "make -j8 USE_CUDA=FALSE USE_MPI=TRUE DEBUG=FALSE;"
    "mpiexec -np 2 ./amr3d.gnu.MPI.ex inputs.3d.FallingSphere max_step=1 amr.n_cell=16 8 8",  
    shell=True, 
    check=True, 
    cwd=working_dir,
    stdout=None if print_output else subprocess.DEVNULL,
    stderr=None if print_output else subprocess.DEVNULL
    )
    print("test_FallingSphere succceed")
    
def test_FlowPastCylinder(working_dir, print_output):
    subprocess.run(
    "make -j8 USE_CUDA=FALSE USE_MPI=TRUE DEBUG=FALSE;"
    "mpiexec -np 2 ./amr3d.gnu.MPI.ex inputs.3d.flow_past_cylinder-x max_step=1",  
    shell=True, 
    check=True, 
    cwd=working_dir,
    stdout=None if print_output else subprocess.DEVNULL,
    stderr=None if print_output else subprocess.DEVNULL
    )
    print("test_FlowPastCylinder succceed")

def test_FlowPastSphere(working_dir, print_output):
    subprocess.run(
    "make -j8 USE_CUDA=FALSE USE_MPI=TRUE DEBUG=FALSE;"
    "mpiexec -np 2 ./amr3d.gnu.MPI.ex inputs.3d.flow_past_sphere max_step=1 amr.n_cell=16 8 8",  
    shell=True, 
    check=True, 
    cwd=working_dir,
    stdout=None if print_output else subprocess.DEVNULL,
    stderr=None if print_output else subprocess.DEVNULL
    )
    print("test_FlowPastSphere succceed")

def test_RayleighTaylor(working_dir, print_output):
    subprocess.run(
    "make -j8 USE_CUDA=FALSE USE_MPI=TRUE DEBUG=FALSE && mpiexec -np 2 ./amr2d.gnu.MPI.ex inputs.2d.rayleightaylor", 
    shell=True, 
    check=True, 
    cwd=working_dir,
    stdout=None if print_output else subprocess.DEVNULL,
    stderr=None if print_output else subprocess.DEVNULL
    )
    print("test_RayleighTaylor succceed")

def test_RayleighTaylor_LS(working_dir, print_output):
    subprocess.run(
    "make -j8 USE_CUDA=FALSE USE_MPI=TRUE DEBUG=FALSE;"
    "mpiexec -np 2 ./amr2d.gnu.MPI.ex inputs.2d.rayleightaylor_rt ", 
    shell=True, 
    check=True, 
    cwd=working_dir,
    stdout=None if print_output else subprocess.DEVNULL,
    stderr=None if print_output else subprocess.DEVNULL
    )
    print("test_RayleighTaylor_LS succceed")



def main():
    # if print_output = fasle，the information of compile and running doesn't display in termianl 
    print_output =False 

    script_dir = os.path.dirname(os.path.abspath(__file__))  
    print("Script Directory:", script_dir)

    # LidDrivenCavity
    working_dir = os.path.join(script_dir, "../Tutorials/LidDrivenCavity")
    print("Test Working Directory:", os.path.abspath(working_dir))
    test_lidDrivenCanvity(working_dir, print_output)

    # LidDrivenCavitySphere
    working_dir = os.path.join(script_dir, "../Tutorials/LidDrivenCavitySphere")
    print("Test Working Directory:", os.path.abspath(working_dir))
    test_LidDrivenCavitySphere(working_dir, print_output)

    # RSV
    working_dir = os.path.join(script_dir, "../Tutorials/RSV")
    print("Test Working Directory:", os.path.abspath(working_dir))
    test_RSV(working_dir, False)
    
    
    # DraftingKissingTumbling  
    working_dir = os.path.join(script_dir, "../Tutorials/DraftingKissingTumbling")
    print("Test Working Directory:", os.path.abspath(working_dir))
    test_DraftingKissingTumbling(working_dir, print_output)

    # FallingSphere 
    working_dir = os.path.join(script_dir, "../Tutorials/FallingSphere")
    print("Test Working Directory:", os.path.abspath(working_dir))
    test_FallingSphere(working_dir, print_output)

    # FlowPastCylinder
    working_dir = os.path.join(script_dir, "../Tutorials/FlowPastCylinder")
    print("Test Working Directory:", os.path.abspath(working_dir))
    test_FlowPastCylinder(working_dir, print_output)

    # FlowPastSphere 
    working_dir = os.path.join(script_dir, "../Tutorials/FlowPastSphere")
    print("Test Working Directory:", os.path.abspath(working_dir))
    test_FlowPastSphere(working_dir, print_output)

    # RayleighTaylor
    working_dir = os.path.join(script_dir, "../Tutorials/RayleighTaylor")
    print("Test Working Directory:", os.path.abspath(working_dir))
    test_RayleighTaylor(working_dir, print_output)    

    # RayleighTaylor_LS
    working_dir = os.path.join(script_dir, "../Tutorials/RayleighTaylor_LS")   
    print("Test Working Directory:", os.path.abspath(working_dir))
    test_RayleighTaylor_LS(working_dir, print_output)    



if __name__== "__main__" :
    main()

