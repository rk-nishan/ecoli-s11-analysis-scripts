"""Build or execute the pinned S11 ANI/PGAP command on an existing Linux runtime."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
FROZEN_RUN = {'release': '2026-06-18.build8602',
 'launcher_sha256': '4f526a2afea4a127503675d68ad402d4de1e3bb57b1bb1001dc48647b8fe641a',
 'container_digest': 'ncbi/pgap@sha256:1f30946971f7285a726f9e19550427dcdb608ce924fe68233d677bf66f2c6dd4',
 'assembly_sha256': '406a7003fa816dd990fae93d2a9148659b35d2b56c21667c2322e3f65478372f'}
INPUT_TEMPLATES = {'input_PGAP_S11_own_use.yaml': 'fasta:\n'
                                '    class: File\n'
                                '    location: S11_contigs.fasta\n'
                                'submol:\n'
                                '    class: File\n'
                                '    location: submol_PGAP_S11.yaml\n'
                                'xpath_fail_initial_asnvalidate: >\n'
                                '    //*[\n'
                                '        ( @severity="ERROR" or @severity="REJECT" )\n'
                                '        and not(contains(@code, '
                                '"GENERIC_MissingPubRequirement")) \n'
                                '        and not(contains(@code, '
                                '"SEQ_DESCR_ChromosomeLocation")) \n'
                                '        and not(contains(@code, "SEQ_DESCR_MissingLineage")) \n'
                                '        and not(contains(@code, "SEQ_DESCR_NoTaxonID")) \n'
                                '        and not(contains(@code, '
                                '"SEQ_DESCR_OrganismIsUndefinedSpecies"))\n'
                                '        and not(contains(@code, '
                                '"SEQ_DESCR_StrainWithEnvironSample"))\n'
                                '        and not(contains(@code, '
                                '"SEQ_DESCR_BacteriaMissingSourceQualifier"))\n'
                                '        and not(contains(@code, '
                                '"SEQ_DESCR_UnwantedCompleteFlag")) \n'
                                '        and not(contains(@code, '
                                '"SEQ_FEAT_BadCharInAuthorLastName")) \n'
                                '        and not(contains(@code, "SEQ_FEAT_ShortIntron")) \n'
                                '        and not(contains(@code, "SEQ_INST_InternalNsInSeqRaw")) \n'
                                '        and not(contains(@code, '
                                '"SEQ_INST_ProteinsHaveGeneralID")) \n'
                                '        and not(contains(@code, "SEQ_PKG_NucProtProblem")) \n'
                                '        and not(contains(@code, '
                                '"SEQ_PKG_ComponentMissingTitle")) \n'
                                '    \n'
                                '    and not(contains(@code, "GENERIC_BadSubmissionAuthorName"))\n'
                                '    and not(contains(@code, "GENERIC_BadFirstName"))\n'
                                '    and not(contains(@code, "GENERIC_BadLastName"))\n'
                                '    ]\n'
                                'xpath_fail_final_asnvalidate: >\n'
                                '    //*[( @severity="ERROR" or @severity="REJECT" )\n'
                                '        and not(contains(@code, '
                                '"GENERIC_MissingPubRequirement")) \n'
                                '        and not(contains(@code, '
                                '"SEQ_DESCR_ChromosomeLocation")) \n'
                                '        and not(contains(@code, "SEQ_DESCR_MissingLineage")) \n'
                                '        and not(contains(@code, "SEQ_DESCR_NoTaxonID")) \n'
                                '        and not(contains(@code, '
                                '"SEQ_DESCR_OrganismIsUndefinedSpecies"))\n'
                                '        and not(contains(@code, '
                                '"SEQ_DESCR_StrainWithEnvironSample"))\n'
                                '        and not(contains(@code, '
                                '"SEQ_DESCR_BacteriaMissingSourceQualifier"))\n'
                                '        and not(contains(@code, '
                                '"SEQ_DESCR_UnwantedCompleteFlag")) \n'
                                '        and not(contains(@code, '
                                '"SEQ_FEAT_BadCharInAuthorLastName")) \n'
                                '        and not(contains(@code, "SEQ_FEAT_ShortIntron")) \n'
                                '        and not(contains(@code, "SEQ_INST_InternalNsInSeqRaw")) \n'
                                '        and not(contains(@code, '
                                '"SEQ_INST_ProteinsHaveGeneralID")) \n'
                                '        and not(contains(@code, '
                                '"SEQ_PKG_ComponentMissingTitle")) \n'
                                '        and not(contains(@code, "SEQ_PKG_NucProtProblem")) \n'
                                '    \n'
                                '    and not(contains(@code, "GENERIC_BadSubmissionAuthorName"))\n'
                                '    and not(contains(@code, "GENERIC_BadFirstName"))\n'
                                '    and not(contains(@code, "GENERIC_BadLastName"))\n'
                                '    and not(contains(@code, "SEQ_INST_BadSeqIdFormat"))\n'
                                '    ]\n',
 'submol_PGAP_S11.yaml': 'organism:\n    genus_species: "Escherichia coli"\n    strain: "S11"\n'}


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--launcher", type=Path, required=True, help="Existing pinned NCBI pgap.py")
    cli.add_argument("--reference-dir", type=Path, required=True, help="Existing PGAP_INPUT_DIR containing the pinned data")
    cli.add_argument("--assembly", type=Path, required=True)
    cli.add_argument("--outdir", type=Path, default=ROOT / "outputs/pgap_rerun")
    cli.add_argument("--docker", default="docker", help="Docker executable or compatible site wrapper")
    cli.add_argument("--execute", action="store_true", help="Run the workflow after the input and runtime checks")
    args = cli.parse_args()
    record = FROZEN_RUN
    for path, expected in [(args.assembly, record["assembly_sha256"]), (args.launcher, record["launcher_sha256"])]:
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Pinned input/source checksum differs: {path.name}")
    target = args.outdir.resolve()
    if target.exists():
        raise FileExistsError("Choose a fresh output directory")
    if args.execute:
        if sys.platform != "linux":
            raise SystemExit("PGAP execution requires an existing compatible Linux/Docker runtime")
        for marker in [".pgap_complete", ".ani_complete"]:
            if not (args.reference_dir / marker).is_file():
                raise ValueError(f"Reference completion marker missing: {marker}")
        subprocess.run([args.docker, "image", "inspect", record["container_digest"]],
                       check=True, stdout=subprocess.DEVNULL)
    target.mkdir(parents=True)
    inputs = target / "input"
    inputs.mkdir()
    for name in ["input_PGAP_S11_own_use.yaml", "submol_PGAP_S11.yaml"]:
        (inputs / name).write_bytes(INPUT_TEMPLATES[name].encode("utf-8"))
    shutil.copyfile(args.assembly, inputs / "S11_contigs.fasta")
    argv = [sys.executable, "-u", str(args.launcher.resolve()), "--use-version", record["release"],
            "--no-self-update", "-n", "--docker", args.docker, "--container-path", record["container_digest"],
            "--cpus", "8", "--memory", "28g", "--taxcheck", "-o", str(target / "results"),
            str(inputs / "input_PGAP_S11_own_use.yaml")]
    plan = {"argv": argv, "assembly_sha256": record["assembly_sha256"],
            "PGAP_INPUT_DIR": str(args.reference_dir.resolve()), "execution_requested": args.execute,
            "release": record["release"], "submission_ready": False,
            "validation_configuration": "Preserved analysis-use author-name and temporary protein-ID exemptions; the embedded YAML preserves the actual author-name and temporary protein-ID validation exemptions",
            "reference_integrity": "Completion markers are checked at execution; prepare the exact release using the official NCBI PGAP instructions"}
    (target / "command_plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    if args.execute:
        temp = target / "tmp"
        temp.mkdir()
        env = dict(os.environ, PGAP_INPUT_DIR=str(args.reference_dir.resolve()), TMPDIR=str(temp))
        with (target / "launcher.log").open("w") as log:
            completed = subprocess.run(argv, cwd=target, env=env, stdout=log, stderr=subprocess.STDOUT)
        (target / "exit_code.json").write_text(json.dumps({"exit_code": completed.returncode}) + "\n")
        raise SystemExit(completed.returncode)
    print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    main()
