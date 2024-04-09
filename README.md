# ReLERNN-pipeline
This repository demonstrates a pipeline to install and run ReLERNN (https://doi.org/10.1093/molbev/msaa038), a deep learning method to produce genome-wide landscapes of recombination. This pipeline is specific for running ReLERNN on the NC State BRC Cluster. 

* Please note that ReLERNN is extrememly sensitive to versions of the software used. We recommend following this installation guide exactly because these versions are confirmed to work with ReLERNN.
* Some of these steps will only need to be done once for installation and some of them will need to be repeated each time you run ReLERNN. We have noted which steps will need to be repeated upon running ReLERNN again. 

## Installation and Example Testing
1. Start an interactive session on node94 for one hour using all the cores available 
   ```
   srun --nodes=1 --ntasks-per-node=36 --time=01:00:00 --partition gpu --nodelist node94 --pty bash -i
   ```
    * Note that node94 is the only updated GPU node that will work with ReLERNN
    * Warning: If you begin the interactive session and "(base)" appears before your username@node94, you have an activated conda environment. This will interfere with the installation. Please run ```conda deactivate``` before proceeding.
    * For testing and running the small example that is provided by ReLERNN, we will use an interactive session on the cluster. For bigger jobs, we submit the job using sbatch. Those instructions will be provided below this section. 
    
2. Load all previously installed modules with correct version numbers
    ```
    module load python/3.10.14
    module load cuda/11.8
    module load cudnn/8.9.7.29_cuda11.x
    ```
    * Note that this step needs to be done every time before you run ReLERNN
      
3. Install additional packages using pip
    ```
    python3 -m pip install --user --upgrade pip
    python3 -m pip install --user wheel
    python3 -m pip install --user --pre --upgrade tensorrt==8.6.1
    python3 -m pip install --user tensorflow[and-cuda]==2.13.0
    ```
4. Manually install the files for Tensor RT specific for Cuda11 and v8.6.1
    ```
    wget https://developer.nvidia.com/downloads/compute/machine-learning/tensorrt/secure/8.6.1/tars/TensorRT-8.6.1.6.Linux.x86_64-gnu.cuda-11.8.tar.gz
    tar -xf TensorRT-8.6.1.6.Linux.x86_64-gnu.cuda-11.8.tar.gz
    ```
5. Test that your environment is ready and working. Make sure Tensor RT works with available GPU

    Set the environment:
    ```
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HOME}/TensorRT-8.6.1.6/lib
    ```
    Test the environment:
    ```
    python3 -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"
    ```
    * Note that these steps needs to be done every time before you run ReLERNN
    * If everything works fine, you'll see an optimization warning that can be ignored and an output message that reads `[PhysicalDevice(name='/physical_device:GPU:0', device_type='GPU')]`. This is confirming that you have access to the GPU. 
      
6. Clone ReLERNN repository in your space on the BRC cluster
      * To do this step, please open another window and login to the cluster, but do not enter into an interactive session. This step must be done on the head/login node. This step is currently set up to install ReLERNN in your home directory. If you want to install in another directory, you will need to change the paths accordingly.

      ```
      git clone https://github.com/kr-colab/ReLERNN.git
      ```
      
      * Now, switch back to the window with the interactive session to run the remainder of the commands
      ```
      cd ReLERNN/
      pip install .
      ```
7. Test the ReLERNN example
      ```
      cd examples/
      ./example_pipeline.sh
      ```
      * Running the example should only take a few minutes. You can watch the progress be printed out. When finished running, type ```exit``` to exit the interactive session
      * The output of the example run will be located in a directory called "example_output". The file containing the final results is called "example.PREDICT.BSCORRECTED.txt"

