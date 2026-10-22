"""Tests of scripts/pipeline_web.py (and arbre.py, qc.py) that need no bioinformatics tool: genome checks, N50, run
history, EasyCGTree names, CheckM2's report and the pairwise summary. Run with: python tests/test_pipeline_inputs.py  (or pytest)."""
import contextlib
import io
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import arbre  # noqa: E402
import pipeline_web as pw  # noqa: E402
import qc  # noqa: E402

ADN = "ACGTACGTGGCC" * 100          # 1 200 bp, 8 G/C sur 12 bases


def ecrire(dossier, nom, contenu):
    chemin = os.path.join(dossier, nom)
    with open(chemin, "w") as fh:
        fh.write(contenu)
    return chemin


def echec(fonction, *args):
    """Lance une fonction qui doit arreter le pipeline ; renvoie le message affiche."""
    sortie = io.StringIO()
    with contextlib.redirect_stdout(sortie):
        try:
            fonction(*args)
        except SystemExit:
            return sortie.getvalue()
    raise AssertionError("the pipeline should have stopped")


def silencieux(fonction, *args):
    with contextlib.redirect_stdout(io.StringIO()):
        return fonction(*args)


def test_statistiques_et_n50():
    with tempfile.TemporaryDirectory() as d:
        # trois contigs de 1 200, 600 et 300 bp, avec des N qui ne comptent pas dans le %GC
        f = ecrire(d, "g.fna", ">c1\n%s\n>c2\n%s\n>c3\n%s\n" % (ADN, ADN[:600], "N" * 300))
        total, gc, acgt, sequences, _, n50 = pw.lire_fasta(f)
        assert (total, sequences, n50) == (2100, 3, 1200)
        assert acgt == 1800 and gc == 1200         # les N ne comptent ni dans le GC ni dans le total A/C/G/T


def test_fichiers_refuses():
    with tempfile.TemporaryDirectory() as d:
        for nom, contenu, message in (("vide.fna", "", "contains no sequence"),
                                      ("entetes.fna", ">a\n>b\n", "contains no sequence"),
                                      ("texte.fna", "hello world\n", "not a FASTA file"),
                                      ("proteine.fna", ">p\nMKVLAAGIVE\n", "not DNA")):
            assert message in echec(pw.lire_fasta, ecrire(d, nom, contenu)), nom


def test_liste_des_genomes():
    with tempfile.TemporaryDirectory() as d:
        for nom in ("A.fna", "B.FNA", "C.fasta", "D.FAS", "E.fna"):
            ecrire(d, nom, ">x\n%s\n" % ADN)
        ecrire(d, "notes.txt", "pas un genome")
        genomes = silencieux(pw.lister_genomes, d, ["qc", "arbre"])
        assert sorted(genomes) == ["A", "B", "C", "D", "E"]      # extension en majuscules comprise, .txt ignore


def test_nombre_minimum():
    with tempfile.TemporaryDirectory() as d:
        for nom in ("A.fna", "B.fna"):
            ecrire(d, nom, ">x\n%s\n" % ADN)
        assert len(silencieux(pw.lister_genomes, d, ["qc", "ani"])) == 2
        assert "needs at least 5 genomes" in echec(pw.lister_genomes, d, ["qc", "arbre"])
        os.remove(os.path.join(d, "B.fna"))
        assert "at least 2 genomes" in echec(pw.lister_genomes, d, ["ani"])


def test_noms_refuses():
    with tempfile.TemporaryDirectory() as d:
        ecrire(d, "A.fna", ">x\n%s\n" % ADN)
        ecrire(d, "Bacillus subtilis (168).fna", ">x\n%s\n" % ADN)
        assert "Bacillus_subtilis_168" in echec(pw.lister_genomes, d, ["ani"])
    with tempfile.TemporaryDirectory() as d:
        ecrire(d, "A.fna", ">x\n%s\n" % ADN)
        ecrire(d, "A.fas", ">x\n%s\n" % ADN)
        assert "same genome name" in echec(pw.lister_genomes, d, ["ani"])


def test_noms_easycgtree():
    assert arbre.nom_easycgtree("GCF_000009045.1_ASM904v1_genomic") == "GCF-000009045.1-ASM904v1-genomic"
    assert arbre.nom_easycgtree("Bacillus-subtilis-168") == "Bacillus-subtilis-168"


def test_qc_et_checkm2():
    with tempfile.TemporaryDirectory() as d:
        genomes, res = os.path.join(d, "g"), os.path.join(d, "results")
        os.makedirs(genomes)
        os.makedirs(os.path.join(res, "checkm2_out_fna"))
        ecrire(genomes, "A.fna", ">x\n%s\n" % ADN)
        ecrire(genomes, "B.fna", ">x\n%s\n" % ADN)
        silencieux(qc.ecrire_qc, genomes, res)                 # sans rapport de CheckM2 : NA
        with open(os.path.join(res, "qc.txt")) as fh:
            assert fh.read().splitlines()[1].split("\t")[-2:] == ["NA", "NA"]
        # colonnes reperees par leur nom, meme dans un autre ordre ; B absent du rapport reste a NA
        ecrire(res, "checkm2_out_fna/quality_report.tsv",
               "Name\tContamination\tCompleteness\tCompleteness_Model_Used\nA\t0.15\t99.99\tNeural Network\n")
        silencieux(qc.ecrire_qc, genomes, res)
        with open(os.path.join(res, "qc.txt")) as fh:
            lignes = [l.split("\t") for l in fh.read().splitlines()]
        assert lignes[1][0] == "A" and lignes[1][2] == "66.7" and lignes[1][-2:] == ["99.99", "0.15"]
        assert lignes[2][0] == "B" and lignes[2][-2:] == ["NA", "NA"]


def test_resultats_gardes_ou_archives():
    with tempfile.TemporaryDirectory() as racine:
        g1, g2 = os.path.join(racine, "g1"), os.path.join(racine, "g2")
        os.makedirs(g1)
        os.makedirs(g2)
        for nom in ("A.fna", "B.fna"):
            ecrire(g1, nom, ">x\n%s\n" % ADN)
        ecrire(g2, "A.fna", ">x\n%s\n" % ADN)
        ecrire(g2, "C.fna", ">x\n%s\n" % ADN[:900])
        res, prot = os.path.join(racine, "results"), os.path.join(racine, "proteins")

        def lancer(dossier):
            genomes = silencieux(pw.lister_genomes, dossier, ["ani"])
            return silencieux(pw.preparer_dossiers, genomes, racine, res, prot)

        entree, fichiers = lancer(g1)
        assert sorted(os.listdir(entree)) == ["A.fna", "B.fna"] and all(os.path.islink(f) for f in fichiers)
        ecrire(res, "ani_resultats.txt", "ancien")
        lancer(g1)                                             # memes genomes : le resultat reste
        assert os.path.exists(os.path.join(res, "ani_resultats.txt"))
        lancer(g2)                                             # autres genomes : il part dans previous_runs/
        assert not os.path.exists(os.path.join(res, "ani_resultats.txt"))
        archives = os.listdir(os.path.join(racine, "previous_runs"))
        assert len(archives) == 1
        assert os.path.exists(os.path.join(racine, "previous_runs", archives[0], "results", "ani_resultats.txt"))


def test_synthese():
    with tempfile.TemporaryDirectory() as res:
        # A et B : meme espece par l'ANI et le dDDH ; A et C : especes differentes, meme genre ;
        # B et C : FastANI n'aligne pas (absente du fichier) ; A et D : ANI >= 95 % sur 30 % du genome seulement
        ecrire(res, "ani_resultats.txt",
               "g/A.fna\tg/B.fna\t98.0\t90\t100\ng/B.fna\tg/A.fna\t97.8\t95\t100\n"
               "g/A.fna\tg/C.fna\t85.0\t60\t100\ng/C.fna\tg/A.fna\t85.2\t62\t100\n"
               "g/A.fna\tg/D.fna\t96.0\t30\t100\ng/D.fna\tg/A.fna\t96.0\t31\t100\n")
        os.makedirs(os.path.join(res, "dddh_out"))
        matrice = ("\tA\tB\tC\tD\nA\t0.0\t0.02\t0.15\t0.03\nB\t0.02\t0.0\t0.16\tNA\n"
                   "C\t0.15\t0.16\t0.0\tNA\nD\t0.03\tNA\tNA\t0.0\n")
        ecrire(res, "dddh_out/dddh_distance.tsv", matrice)
        ecrire(res, "dddh_out/dddh_similarite.tsv", matrice)
        os.makedirs(os.path.join(res, "pocp_out"))
        ecrire(res, "pocp_out/pocp_matrice.tsv",
               "\tA\tB\tC\tD\nA\t100\t95\t60\t90\nB\t95\t100\t40\t90\nC\t60\t40\t100\t50\nD\t90\t90\t50\t100\n")
        silencieux(pw.ecrire_synthese, res)
        with open(os.path.join(res, "synthese_paires.tsv")) as fh:
            lignes = [l.rstrip("\n").split("\t") for l in fh]
        col = {nom: i for i, nom in enumerate(lignes[0])}
        paires = {(l[0], l[1]): l for l in lignes[1:]}
        assert len(paires) == 6
        ab, ac, ad, bc = paires[("A", "B")], paires[("A", "C")], paires[("A", "D")], paires[("B", "C")]
        assert ab[col["ani"]] == "97.90" and ab[col["ani_aligned_fraction"]] == "0.90"
        assert ab[col["same_species"]] == "same" and ab[col["same_genus"]] == "same"
        assert ac[col["same_species"]] == "different" and ac[col["same_genus"]] == "same"
        assert ad[col["same_species"]] == "same (low aligned fraction)"
        assert bc[col["ani"]] == "NA" and bc[col["same_species"]] == "different"
        assert bc[col["same_genus"]] == "different" and bc[col["aai"]] == "NA"


if __name__ == "__main__":
    for nom, fonction in sorted(globals().items()):
        if nom.startswith("test_"):
            fonction()
            print("ok  ", nom)
