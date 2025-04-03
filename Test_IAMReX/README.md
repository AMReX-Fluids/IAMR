Test_IAMReX is a custom set of Python scripts designed to check the completeness of the code and verify the correctness of the results it produces. We provided some test cases, including lidDrivenCanvity, RSV, DraftingKissingTumbling, FlowPastCylinder, FlowPastSphere and RayleighTaylor. These tests cover various aspects and involve multiple modules. Additionally, you can follow the provided examples to create tests for yourself.

Meanwhile, with the help of GitHub Actions, we have set up a CI pipeline to automatically run and validate these tests, ensuring continuous integration and testing of the codebase.

# How to use 

test_IAMReX.py is all we need and you can just run it to test some given cases. With the help of GitHub Actions, this process can be automated. You can find the CI task in .github/workflows/test_IAMReX.yml.

```bash
python3 test_IAMReX.py


# Here are the running results, which indicate that the case has been built and run successfully.
Script Directory: xxx/IAMReX/Test_IAMReX
Test Working Directory: xxx/IAMReX/Tutorials/RayleighTaylor_LS
test_RayleighTaylor_LS succceed
```

If you want to add a new test, here are only two thing you need to do. 

1. Define the test function as follow. A test is defined here to check the compilation and execution of the Rayleigh-Taylor-related code. The working_diry is a relative path to test_IAMReX.py. It is used to locate and execute the necessary files for testing.

    ```python
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
    ```
2. Add it to the main() .
    ```python
    def main():
        # RayleighTaylor
        working_dir = os.path.join(script_dir, "../Tutorials/RayleighTaylor")
        print("Test Working Directory:", os.path.abspath(working_dir))
        test_RayleighTaylor(working_dir, print_output)   
    ```



NOTE: 

For github action, mpi is only supported 2 cores, GPU is not supported.


