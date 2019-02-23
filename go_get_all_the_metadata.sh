#go_get_all_the_metadata.sh

cd /vol/sci/astro/cosmo/nas2/NAS2/sroca/Simulations/NIHAO

for d in g*//; do
    echo "$d"
	cd ~/quasarscan/galaxy_catalogs_nihao
	mkdir "$d"
    cd /vol/sci/astro/cosmo/nas2/NAS2/sroca/Simulations/NIHAO
    cp "$d"/analysis/Mvir_fbar_R200.out ~/quasarscan/galaxy_catalogs_nihao/"$d"
done