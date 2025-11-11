function out = egi1010(in,search)
% egi1020() - > Returns a dictionary mapping 1010 labels to EGI channels
% egi1020("C3") -> Returns the channel corresponding to C3
% egi1020([34 44]) -> returns the labels corrsponding to these 2 channels
% egi("O",@startsWith)  -> Returns all channels with a label that starts with O
%
arguments
    in =[]
    search = []  % Use @startsWith or @contains to search for groups of labels.
end

% Data from https://www.egi.com/images/HydroCelGSN_10-10.pdf
T = cell2table({"C3"	59	1.029910
    "AF3"	34	0.509761
    "C4"	183	0.623156
    "C1"	44	0.519046
    "PO8"	161	0.525352
    "F3"	36	0.685314
    "AF4"	12	0.526048
    "F4"	224	0.590697
    "Afz"	20	0.528535
    "F7"	47	0.537374
    "TP8"	179	0.548095
    "F8"	2	1.343284
    "FC3"	42	0.561542
    "FP1"	37	1.053985
    "CP3"	66	0.589000
    "FP2"	18	0.811629
    "P6"	162	0.638474
    "FPZ"	26	1.218066
    "PO3"	109	0.677353
    "Fz"	21	0.845529
    "C2"	185	0.704042
    "O1"	116	1.257897
    "FC1"	24	0.704503
    "O2"	150	0.959035
    "PO4"	140	0.759366
    "P3"	87	0.467009
    "Oz"	126	0.766504
    "P4"	153	0.613058
    "Cp2"	143	0.768954
    "T3"	69	0.565707
    "T7"	69	0.565707
    "FC2"	207	0.771819
    "T4"	202	0.347099
    "T8"	202	0.347099
    "CP1"	79	0.800214
    "T5"	96	0.766165
    "P7"	96	0.766165
    "TP9"	94	0.801732
    "T6"	170	0.394401
    "P8"	170	0.394401
    "F1"	29	0.865540
    "Pz"	101	0.126808
    "Fcz"	15	0.916188
    "Poz"	119	0.158625
    "TP10"	190	0.917598
    "F2"	5	0.160534
    "F10"	226	0.922800
    "FC5"	49	0.227681
    "P2"	142	1.012270
    "Ft10"	219	0.242712
    "F5"	48	1.012833
    "C6"	194	0.253742
    "P9"	106	1.048151
    "Ft9"	67	0.259169
    "FC4"	206	1.079352
    "F6"	222	0.298570
    "Cp5"	76	1.084715
    "FT8"	211	0.316283
    "FC6"	213	1.120079
    "AF8"	10	0.321524
    "Afz"	27	1.137608
    "CpZ"	81	0.326043
    "Po7"	97	1.193410
    "CP6"	172	0.344188
    "AF7"	46	1.211161
    "C5"	64	0.354355
    "Afz"	26	1.338032
    "CP4"	164	0.382313
    "Tp7"	84	1.344355
    "P10"	169	0.432754
    "FT7"	62	1.365872
    "F9"	252	0.443114
    "T9"	68	3.302035
    "P1"	88	0.464150
    "T10"	210	3.434444
    "P5"	86	0.500231},'VariableNames',{'Label','Channel','Distance'});


D = dictionary(T.Label,T.Channel);
if ~isempty(search)
    out = T.Channel(search(T.Label,in));
elseif isempty(in)
    % Return the full dictionary
    out =D;
elseif isnumeric(in) && all(ismember(in,T.Channel))
    % Lookup the corresponding 1020 Label
    out = T.Label(ismember(in,T.Channel));
elseif (isstring(in) || ischar(in) ) && all(ismember(string(in),T.Label))
    % Reverse lookup - return the channel
    out = T.Channel(ismember(upper(T.Label),upper(string(in))));
else
    in
    error('egi1020 needs either an empty input (returns the dictionary), a channel number (Returns the 1020 name) or a 1020 label (Returns the channel number)')
end


