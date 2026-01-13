# c-morpho-butterfly

Hamid Naderi Yeganeh is an artist which draws pictures with mathematical equations: [Wikipedia](https://en.wikipedia.org/wiki/Hamid_Naderi_Yeganeh), [X](https://x.com/naderi_yeganeh), [Instagram](https://www.instagram.com/hamidnaderiyeganeh/), [YouTube](https://www.youtube.com/@naderiyeganeh).

Morpho buttterfly is one of his artworks ([X](https://www.instagram.com/p/DQmRlpliAVM/), [Instagram](https://x.com/naderi_yeganeh/status/1754107434224804305/)):
![original butterfly and equations](pics/original-morpho-butterfly-and-equations.jpg)

The present repository reproduces this work with C.

The program uses OpenMP to parallelize calculations and [STB library](https://github.com/nothings/stb) to create the images.

The program:  
\- generates butterfly picture (by changing `NB_RUNS`, it can be done several times, to evaluate durations), and  
\- generates heatmaps of intermediate calculations.

See images below.

Any comment? Open an [issue](https://github.com/occisn/c-morpho-butterfly/issues), or start a discussion [here](https://github.com/occisn/c-morpho-butterfly/discussions) or [at profile level](https://github.com/occisn/occisn/discussions).

([the same in Common Lisp](https://github.com/occisn/cl-morpho-butterfly))

Generated butterfly:  
![output butterfly](pics/output-butterfly.jpg)

C heatmap:  
![C heatmap](pics/C-heatmap.jpg)

E heatmap:  
![E heatmap](pics/E-heatmap.jpg)

L heatmap:  
![L heatmap](pics/L-heatmap.jpg)

W heatmap:  
![L heatmap](pics/W-heatmap.jpg)

A0 heatmap:  
![A0 heatmap](pics/A0-heatmap.jpg)

A1 heatmap:  
![A1 heatmap](pics/A1-heatmap.jpg)

K0 heatmap:  
![K0 heatmap](pics/K0-heatmap.jpg)

K1 heatmap:  
![K1 heatmap](pics/K1-heatmap.jpg)

K2 heatmap:  
![K2 heatmap](pics/K2-heatmap.jpg)

H0 heatmap:  
![H0 heatmap](pics/H0-heatmap.jpg)

H1 heatmap:  
![H1 heatmap](pics/H1-heatmap.jpg)

H2 heatmap:  
![H2 heatmap](pics/H2-heatmap.jpg)

(end of README)
