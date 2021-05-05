#!/bin/bash -l


exec_AugmentPlotfile='AugmentPlotfile3d.gnu.MPI.ex'
exec_DeriveSpectrum='AmrDeriveSpectrum3d.gnu.MPI.ex'

# Check Needed Files are in Current Directory
reqfiles_ascii='yt_post_scrape_tseries.py yt_post_scrape_common.py yt_post_scrape_slice.py yt_post_scrape_common.py spectra.py derivespect-inputs '
reqfiles=$reqfiles_ascii $exec_DeriveSpectrum $exec_AugmentPlotfile
for f in $reqfiles; do
    if [ ! -f $f ]; then
        echo "File $f needed by script not found. Exiting..."
        exit 1
    fi
done



# Check Required Arguments Provided
if [ $# -ne 2 ]; then
    echo "Must provide exactly one argument; the name of the folder to be processed, and a flag if you want to augment plotfiles (very time consuming ...)."
    exit
else
    echo "Processing $1..."
fi


nprocs=4
folder=$1
tmpfname='tmp_plotfile'
run="mpirun -n $nprocs"
pyrun="$run python"


if [ $2 -gt 0 ];then
  # Augment Plotfile
  echo "Augmenting Plotfiles..."

  echo "TOTO"
  echo "$exec_AugmentPlotfile"
  echo "$folder"

  cp $exec_AugmentPlotfile $folder

  for f in `ls -1 $folder | grep "plt_[0-9]*[0-9]" | sed '/csv/d' | sed '/pickle/d'`; do
      echo "Processing ${f}..."
      cd $folder
      echo -e "infile = ${f}\noutfile = ${f}\nadd_vorticity = 1\nadd_divergence = 1" > tmp_inputs
      $run ./$exec_AugmentPlotfile tmp_inputs
      rm -r *.old.*
      cd ..
  done
  rm $folder/$exec_AugmentPlotfile
  rm $folder/tmp_inputs
else
  echo "Not Augmenting Plotfiles..."
fi



# Time Series
echo "Computing time series for $folder."
$pyrun yt_post_scrape_tseries.py -m $folder $folder


# Final Slice Processing
last=`ls -1 $folder | grep plt | sed '/csv/d' | sed '/pickle/d' | tail -n 1`
echo -e "Last plotfile:\t$last"

echo "Processing $last..."
mv $folder/$last $tmpfname
$run ./$exec_DeriveSpectrum derivespect-inputs
paste -d ' ' $tmpfname/x_vel*_spectrum.dat $tmpfname/y_vel*_spectrum.dat $tmpfname/z_vel*_spectrum.dat > $tmpfname/vel_spectrum.dat
paste -d ' ' $tmpfname/x_vort*_spectrum.dat $tmpfname/y_vort*_spectrum.dat $tmpfname/z_vort*_spectrum.dat > $tmpfname/vort_spectrum.dat
mv $tmpfname $folder/$last
$pyrun yt_post_scrape_slice.py -d $folder $last

