"""
These rules will only run during validation mode

Takes:
    - list of genome locations
    - table specifying how many reads / fold change to
    simulate for each genome
Output:
    - table of the sensitivity / specificity of the pipline
    - plots of this data
"""

rule simulate_reads:
    input:
        "simulation_input.txt"
    output:
        "00_simulated_reads/sim_R1.fastq",
        "00_simulated_reads/sim_R2.fastq",





