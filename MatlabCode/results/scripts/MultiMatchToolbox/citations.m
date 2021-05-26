%citations.m
%prints out a list of the citations that should be referenced in any
%published work.
%references

s1=('**CITATIONS**\nPlease cite the following references in any publications:\nJarodzka, H.,');
s2=(' Holmqvist, K., & Nyström, M. (2010). A vector-based,\nmultidimensional');
s3=(' scanpath similarity measure. Proceedings of the 2010\nSymposium on Eye-Tracking');
s4=(' Research & Applications, 211-218.\n\nDewhurst, R., Nyström, M., Jarodzka, H., Foulsham, T.,');
s5=(' Johansson, R., & Holmqvist, K. (2012).\nIt depends on how you look at it: Scanpath');
s6=(' comparison in multiple dimensions with Multimatch,\na vector-based approach. Behavior Research');
s7=(' Methods, Springer, 1-22\n\nsee http://wiki.humlab.lu.se/dokuwiki/doku.php?id=public:useful_links&s%%5B%%5D=multimatch\n\n');
refs = strcat(s1,s2,s3,s4,s5,s6,s7);
fprintf(refs);