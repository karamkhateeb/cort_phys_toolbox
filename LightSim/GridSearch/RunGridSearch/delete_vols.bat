echo off
for %%a in (1mm 1p5mm 2mm p5mm) do (
	cd %%a
	for /D %%b in (*) do (
		if exist %%b/*.mc2 if exist %%b/*.nii rm %%b/*.nii)
	cd ..)