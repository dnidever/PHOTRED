\rm /net/home/dln5q/photred/sav/*
echo "Making PHOTRED IDL .sav files" 
idl < mk_photred_sav.batch
idl < mk_photred_rename_sav.batch
idl < mk_photred_split_sav.batch
idl < mk_photred_wcs_sav.batch
idl < mk_photred_daophot_sav.batch
idl < mk_photred_match_sav.batch
idl < mk_photred_allframe_sav.batch
idl < mk_photred_apcor_sav.batch
idl < mk_photred_calib_sav.batch
idl < mk_photred_astrom_sav.batch
idl < mk_photred_combine_sav.batch
idl < mk_photred_deredden_sav.batch
idl < mk_photred_save_sav.batch
idl < mk_photred_html_sav.batch
idl < mk_photred_summary_sav.batch
idl < mk_stdred_sav.batch
idl < mk_stdred_rename_sav.batch
idl < mk_stdred_aperphot_sav.batch
idl < mk_stdred_daogrow_sav.batch
idl < mk_stdred_astrom_sav.batch
idl < mk_stdred_matchcat_sav.batch
idl < mk_stdred_combinecat_sav.batch
idl < mk_stdred_fitdata_sav.batch
idl < mk_stdred_summary_sav.batch
idl < mk_stdred_copy_sav.batch

echo "Making PHOTRED IDL .sav TAR file"
#\rm /net/home/dln5q/photred/photred_idlsav.tar
#cp /net/home/dln5q/idl/factorsum.dat /net/home/dln5q/photred/sav/
#cd sav
#tar -cvf photred_idlsav.tar *.sav factorsum.dat
#\cp photred_idlsav.tar ../
#cd ../

#echo "Copying the TAR file to public_html"
#\cp photred_idlsav.tar /net/home/dln5q/public_html/research/photred/
