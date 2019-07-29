#move_pleiades_files_over.sh

VELANUM=$1

cd /lou/la1/dceverin/VELA_v2.0/
cd VELA$VELANUM
echo "entering directory VELA$VELANUM"

for d in 500 400 330 250 200; do
	echo "copying 10MpcBox_csf512_a0.$d.d..."
    shiftc 10MpcBox_csf512_a0.$d.d /nobackupp2/cstrawn/mydir/VELA${VELANUM}v2
	echo "copying PMcrda0.$d.DAT..."
    shiftc PMcrda0.$d.DAT /nobackupp2/cstrawn/mydir/VELA${VELANUM}v2
	echo "copying PMcrs0a0.$d.DAT..."
    shiftc PMcrs0a0.$d.DAT /nobackupp2/cstrawn/mydir/VELA${VELANUM}v2
	echo "copying pta0.$d.DAT..."
    shiftc pta0.$d.dat /nobackupp2/cstrawn/mydir/VELA${VELANUM}v2
	echo "copying stars_a0.$d.dat..."
    shiftc stars_a0.$d.dat /nobackupp2/cstrawn/mydir/VELA${VELANUM}v2
done