foreach file ( `ls *.c *.C`)
	set ofile = `echo $file | sed 's/\.C/\.o/' | sed 's/\.c/\.o/'`
	set hfile = `echo $file | sed 's/\.C/\.h/' | sed 's/\.c/\.h/'`

	echo $ofile

	if( ! -e $ofile ) then
		mv $file garbage
		if( -e $hfile ) then
			mv $hfile garbage
		endif
	endif
end
