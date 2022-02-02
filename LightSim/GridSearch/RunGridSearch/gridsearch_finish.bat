for %%a in (p5mm 1mm 1p5mm 2mm) do (
	cd %%a
	for /D %%b in (*) do (
		if not exist %%b/*.mc2 "C:\Program Files\MCXStudio\MCXSuite\mcx\bin\mcx.exe" --session "grid" --input %%b/specs.json --root %%b --outputformat mc2 --gpu 1000000 --autopilot 1 --photon 1000000 --normalize 1 --save2pt 1 --reflect 1 --savedet 0 --unitinmm 0.01 --seed 1648335518 --saveseed 0 --momentum 0 --skipradius -2 --array 0 --dumpmask 1 --repeat 1 --maxdetphoton 1000000
		if not exist %%b/*.mc2 "C:\Program Files\MCXStudio\MCXSuite\mcx\bin\mcx.exe" --input %%b/specs.json --root %%b --outputformat mc2 --gpu 1000000 --autopilot 1 --photon 1000000 --normalize 1 --save2pt 1 --reflect 1 --savedet 0 --unitinmm 0.01 --seed 1648335518 --saveseed 0 --momentum 0 --skipradius -2 --array 0 --dumpmask 0 --repeat 1 --maxdetphoton 1000000)
	cd ..)