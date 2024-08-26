echo "Downloading requested GERP data..."

wget -cO https://ftp.ensembl.org/pub/release-112/compara/conservation_scores/91_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
> ../../data/raw/gerp_91mammals_conservation_scores.homo_sapiens.GRCh38.bw

echo "Done."
