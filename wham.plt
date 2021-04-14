# Set the output to a png file
set terminal png size 500,500
# The file we'll write to
set output 'wham_free_eng.png'
# The graphic title
set title 'Free energy vs. Distance(A)'
set xlabel "Distance(A)"
set ylabel "Free Energy (kcal/mol)"
#plot the graphics
plot 'wham.dat' using 1:2 lt rgb "blue" title 'wham' with linespoints, '

