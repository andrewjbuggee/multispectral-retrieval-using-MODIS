% --- My Saved Colors ---

% A hodgepodge of colors I've saved over time because I found tehm
% appealing

% INPUTS
%   n - the number of colors you'd like
%   'random' - the string input tells the function that you want either
%   some precise color, or a random assortment of colors. The input
%   possibilities are:
%       (a) 'random' - function will output a random assortment of n colors
%       (b) 'fixed' - function will output the color codes associated with
%       the vector n according to how they have been saved in the matrix
%       below

function C = mySavedColors(n, randomORfixed)


savedColors = [8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01;...        % (1) a nice bright pink
    4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01;...                   % (2) sea foam green
    0.939001561999887   0.587044704531417   0.230488160211558;...                                   % (3) Kraft Mac and Cheese Orange
    7.065829534842050e-02     9.583877169541491e-01     6.998897395776505e-01;...                   % (4) a bright aquamarine
    0.827450980392157, 0.368627450980392, 0.376470588235294;...                                     % (5) A pleasing salmon red
    3.563645953575846e-01     4.380836048512262e-01     5.147715889386915e-01;...                   % (6) a pleasing blueish gray
    4.672877744560139e-01     8.367186588778170e-01     2.206248068838936e-01;...                   % (7) lime green
    8.259258869082566e-01     6.171862413652410e-01     1.938744022737314e-01;...                   % (8) matte harvest gold
    0.350727103576883   0.622475086001227   0.470923348517591;...                                   % (9) Matte Irish green
    0.550156342898422   0.301246330279491   0.194764289567049;...                                   % (10) UPS brown
    0.80,0.79,0.85;...                                                                              % (11) A pale grey
    0.96,0.42,0.65;...                                                                              % (12) Bubble gum pink
    0.49,0.92,0.04;...                                                                              % (13) Neon green
    0.14,0.96,0.93;                                                                                 % (14) A bright electric blue
    6.017416069441683e-02     4.714776469952955e-01     3.110165739381916e-01;...                   % (15) forest green
    0.471937460337157,0.064164021255792,0.581686629640886;...                                       % (16) Barney purple
    0.264703182911644,0.280822855960028,0.819417369495718;...                                       % (17) a purple-ish blue
    0.087146608906061,0.275650182414309,0.919165103980364;...                                       % (18) a mellow blue
    0.242602763746494,0.165654616103410,0.589586440293093;...                                       % (19) darker purple
    0.097641770151945,0.157654610063639,0.570247296307671;...                                       % (20) a dark blue
    8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01;...                   % (1) a nice bright pink
    4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01;...                   % (2) sea foam green
    0.939001561999887   0.587044704531417   0.230488160211558;...                                   % (3) Kraft Mac and Cheese Orange
    7.065829534842050e-02     9.583877169541491e-01     6.998897395776505e-01;...                   % (4) a bright aquamarine
    0.827450980392157, 0.368627450980392, 0.376470588235294;...                                     % (5) A pleasing salmon red
    3.563645953575846e-01     4.380836048512262e-01     5.147715889386915e-01;...                   % (6) a pleasing blueish gray
    4.672877744560139e-01     8.367186588778170e-01     2.206248068838936e-01;...                   % (7) lime green
    8.259258869082566e-01     6.171862413652410e-01     1.938744022737314e-01;...                   % (8) matte harvest gold
    0.350727103576883   0.622475086001227   0.470923348517591;...                                   % (9) Matte Irish green
    0.550156342898422   0.301246330279491   0.194764289567049;...                                   % (10) UPS brown
    0.80,0.79,0.85;...                                                                              % (11) A pale grey
    0.96,0.42,0.65;...                                                                              % (12) Bubble gum pink
    0.49,0.92,0.04;...                                                                              % (13) Neon green
    0.14,0.96,0.93;                                                                                 % (14) A bright electric blue
    6.017416069441683e-02     4.714776469952955e-01     3.110165739381916e-01;...                   % (15) forest green
    0.471937460337157,0.064164021255792,0.581686629640886;...                                       % (16) Barney purple
    0.264703182911644,0.280822855960028,0.819417369495718;...                                       % (17) a purple-ish blue
    0.087146608906061,0.275650182414309,0.919165103980364;...                                       % (18) a mellow blue
    0.242602763746494,0.165654616103410,0.589586440293093;...                                       % (19) darker purple
    0.097641770151945,0.157654610063639,0.570247296307671;...                                       % (20) a dark blue
    8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01;...                   % (1) a nice bright pink
    4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01;...                   % (2) sea foam green
    0.939001561999887   0.587044704531417   0.230488160211558;...                                   % (3) Kraft Mac and Cheese Orange
    7.065829534842050e-02     9.583877169541491e-01     6.998897395776505e-01;...                   % (4) a bright aquamarine
    0.827450980392157, 0.368627450980392, 0.376470588235294;...                                     % (5) A pleasing salmon red
    3.563645953575846e-01     4.380836048512262e-01     5.147715889386915e-01;...                   % (6) a pleasing blueish gray
    4.672877744560139e-01     8.367186588778170e-01     2.206248068838936e-01;...                   % (7) lime green
    8.259258869082566e-01     6.171862413652410e-01     1.938744022737314e-01;...                   % (8) matte harvest gold
    0.350727103576883   0.622475086001227   0.470923348517591;...                                   % (9) Matte Irish green
    0.550156342898422   0.301246330279491   0.194764289567049;...                                   % (10) UPS brown
    0.80,0.79,0.85;...                                                                              % (11) A pale grey
    0.96,0.42,0.65;...                                                                              % (12) Bubble gum pink
    0.49,0.92,0.04;...                                                                              % (13) Neon green
    0.14,0.96,0.93;                                                                                 % (14) A bright electric blue
    6.017416069441683e-02     4.714776469952955e-01     3.110165739381916e-01;...                   % (15) forest green
    0.471937460337157,0.064164021255792,0.581686629640886;...                                       % (16) Barney purple
    0.264703182911644,0.280822855960028,0.819417369495718;...                                       % (17) a purple-ish blue
    0.087146608906061,0.275650182414309,0.919165103980364;...                                       % (18) a mellow blue
    0.242602763746494,0.165654616103410,0.589586440293093;...                                       % (19) darker purple
    0.097641770151945,0.157654610063639,0.570247296307671;...                                       % (20) a dark blue
    [255,147,215]./255;...                                                                          % (61) colorblind friendly #1
    [114,226,115]./255;...                                                                          % (62) colorblind friendly #2
    [99,0,25]./255;...                                                                              % (63) colorblind friendly #3
    [55,57,169]./255;...                                                                            % (64) colorblind friendly #4
    [238,172,47]./255;...                                                                           % (65) colorblind friendly #5
    [199,171,255]./255;...                                                                          % (66) colorblind friendly #6
    [75,114,0]./255;...                                                                             % (67) colorblind friendly #7
    [255,118,180]./255;...                                                                          % (68) colorblind friendly #8 
    [167,44,1]./255;...                                                                             % (69) colorblind friendly #9
    [0,219,233]./255;...                                                                            % (70) colorblind friendly #10
    [145,127,197]./255;...                                                                          % (71) colorblind friendly #11
    [255,105,101]./255;...                                                                          % (72) colorblind friendly #12
    [217,39,97]./255;...                                                                            % (73) colorblind friendly #13
     ];


if strcmp(randomORfixed, 'random')==true

    indices2use = randperm(size(savedColors,1), n);          % random set of n unique integers between the values 1 and the length of my saved colors.
    C = savedColors(indices2use,:);

elseif strcmp(randomORfixed, 'fixed')==true
    
    C = savedColors(n,:);


else

    error([newline,'I dont know if you want a random assortment of colors or some precise colors', newline])

end





end