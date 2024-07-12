% Mb=7;
% [incli]=deltam(Mb);
% pendiente=incli;

Mb=7;
while Mb<101
    [incli]=deltam(Mb);
    pendiente=incli;
    Mb=Mb+2;
end
