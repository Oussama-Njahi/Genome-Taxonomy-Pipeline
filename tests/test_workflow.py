"""Tests of the Snakemake workflow (workflow/Snakefile) as pipeline_web.py runs it, in dry-run mode (snakemake -n): which
jobs each selection of steps runs. They need Snakemake but no bioinformatics tool. Run with:
python tests/test_workflow.py  (or pytest)."""
import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import pipeline_web as pw  # noqa: E402

ADN = "ACGTACGTGGCC" * 100
NOMS = ["A", "B", "C", "D", "E"]


def taches(racine, steps):
    """Prepare a run of the genomes NOMS in racine as pipeline_web.py does, and return the steps it plans and the
    number of jobs of each rule that Snakemake would run."""
    genomes_dir, res, prot = os.path.join(racine, "g"), os.path.join(racine, "results"), os.path.join(racine, "proteins")
    os.makedirs(genomes_dir, exist_ok=True)
    for nom in NOMS:
        chemin = os.path.join(genomes_dir, nom + ".fna")
        if not os.path.exists(chemin):        # ecrit une seule fois : un genome plus recent qu'un resultat le perime
            with open(chemin, "w") as fh:
                fh.write(">x\n%s\n" % ADN)
    selection = [k for k in pw.STEPS_KEYS if k in steps]
    with contextlib.redirect_stdout(io.StringIO()):
        genomes = pw.lister_genomes(genomes_dir, selection)
        entree, _ = pw.preparer_dossiers(genomes, racine, res, prot)
    prevues = pw.etapes_prevues(selection, list(genomes), prot)
    configfile = os.path.join(res, "snakemake_config.json")
    with open(configfile, "w") as fh:
        json.dump({"genomes": entree, "results": res, "proteins": prot, "steps": prevues,
                   "job_id": "test", "outgroup": None}, fh)
    # sans --sdm conda : l'environnement de CheckM2 n'est pas necessaire pour un essai a blanc
    cmd = [c for c in pw.commande_snakemake(prevues, selection, racine, configfile, 4)
           if c not in ("--sdm", "conda", "--quiet", "rules", "host", "progress")]
    sortie = subprocess.run(cmd + ["-n", "--quiet", "rules"], capture_output=True, text=True)
    assert sortie.returncode == 0, sortie.stdout + sortie.stderr
    lignes = sortie.stdout.split("Job stats:")[1].splitlines()[3:]
    compte = {}
    for ligne in lignes:
        mots = ligne.split()
        if len(mots) != 2 or mots[0] == "total":
            break
        compte[mots[0]] = int(mots[1])
    return prevues, compte


def test_toutes_les_etapes():
    with tempfile.TemporaryDirectory() as d:
        prevues, compte = taches(d, pw.STEPS_KEYS)
        assert prevues == pw.STEPS_KEYS
        assert compte["prodigal"] == compte["ezaai_extract"] == len(NOMS)     # une tache par genome
        for regle in ("qc", "checkm2", "proteines", "ani", "dddh", "aai", "pocp", "arbre", "figures"):
            assert compte[regle] == 1, regle


def test_pocp_lance_prodigal_si_les_proteines_manquent():
    with tempfile.TemporaryDirectory() as d:
        prevues, compte = taches(d, ["pocp"])
        assert prevues == ["proteines", "pocp"] and compte["prodigal"] == len(NOMS)
        for nom in NOMS:                                  # proteines deja la (run precedent) : POCP seul
            with open(os.path.join(d, "proteins", nom + ".faa"), "w") as fh:
                fh.write(">p\nMKV\n")
        prevues, compte = taches(d, ["pocp"])
        assert prevues == ["pocp"] and compte == {"pocp": 1}


def test_checkm2_refait_qc():
    with tempfile.TemporaryDirectory() as d:
        prevues, compte = taches(d, ["checkm2"])
        assert prevues == ["qc", "checkm2"] and compte == {"checkm2": 1, "qc": 1}


def test_figures_seules():
    with tempfile.TemporaryDirectory() as d:
        prevues, compte = taches(d, ["figures"])
        assert prevues == ["figures"] and compte == {"figures": 1}


if __name__ == "__main__":
    for nom, fonction in sorted(globals().items()):
        if nom.startswith("test_"):
            fonction()
            print("ok  ", nom)
