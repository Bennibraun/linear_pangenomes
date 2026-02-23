from pathlib import Path

configfile: "config/config.yaml"

SCRIPTS = config["scripts"]
MARKER_DIR = Path(config.get("marker_dir", "results/markers"))

shell.executable("bash")

rule all:
    input:
        MARKER_DIR / "find_data_sra.done",
        MARKER_DIR / "assemblies.done",
        MARKER_DIR / "call_svs.done",
        MARKER_DIR / "make_cactus_graph.done",
        MARKER_DIR / "align_wgs.done"

rule find_data_sra:
    output:
        MARKER_DIR / "find_data_sra.done"
    conda:
        "envs/base.yaml"
    params:
        script=SCRIPTS["find_data_sra"]
    shell:
        r"""
        set -euo pipefail
        python {params.script}
        mkdir -p {MARKER_DIR}
        touch {output}
        """

rule assemblies:
    input:
        fastq=config["paths"]["assembly_fastq"]
    output:
        MARKER_DIR / "assemblies.done"
    conda:
        "envs/base.yaml"
    params:
        script=SCRIPTS["assemblies"]
    shell:
        r"""
        set -euo pipefail
        bash {params.script} {input.fastq}
        mkdir -p {MARKER_DIR}
        touch {output}
        """

rule call_svs:
    output:
        MARKER_DIR / "call_svs.done"
    conda:
        "envs/base.yaml"
    params:
        script=SCRIPTS["call_svs"]
    shell:
        r"""
        set -euo pipefail
        bash {params.script}
        mkdir -p {MARKER_DIR}
        touch {output}
        """

rule make_cactus_graph:
    output:
        MARKER_DIR / "make_cactus_graph.done"
    conda:
        "envs/base.yaml"
    params:
        script=SCRIPTS["make_cactus_graph"]
    shell:
        r"""
        set -euo pipefail
        bash {params.script}
        mkdir -p {MARKER_DIR}
        touch {output}
        """

rule align_wgs:
    output:
        MARKER_DIR / "align_wgs.done"
    conda:
        "envs/base.yaml"
    params:
        script=SCRIPTS["align_wgs"]
    shell:
        r"""
        set -euo pipefail
        bash {params.script}
        mkdir -p {MARKER_DIR}
        touch {output}
        """
