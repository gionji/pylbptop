
import numpy as np
import dec2bin

#  MAPPING = GETMAPPING(SAMPLES,MAPPINGTYPE) returns a
#  structure containing a mapping table for
#  LBP codes in a neighbourhood of SAMPLES sampling
#  points. Possible values for MAPPINGTYPE are
#       'u2'   for uniform LBP
#       'ri'   for rotation-invariant LBP
#       'riu2' for uniform rotation-invariant LBP.

## get mapping

def sum( n):
    n = str(bin(n))

    one_count = 0
    for i in n:
        if i == "1":
            one_count+=1
    return one_count


def get_mapping(samples, mapping_type):
    table = range(0, 2**samples-1)
    newMax  = 0; #number of patterns in the resulting LBP code
    index   = 0;

    if mappingtype =='u2':
        newMax = samples*(samples-1) + 3;
        for i in range(0, 2**samples-1):
            i_bin = dec2bin(i,samples); # pip install dec 2 bin

            # circularly rotate left number of 1->0 and in binary string
            # x is equal to the number of 1-bits ifmain XOR(x,Rotate left(x))
            j_bin = np.roll(i_bin,-1);
            numt = sum(i_bin ^ j_bin);

            if numt <= 2:
                table(i+1) = index;
                index = index + 1;
            else:
                table(i+1) = newMax - 1;

        if mappingtype == 'ri': #Rotation invariant
            tmpMap = zeros(2**samples,1) - 1;

            for i in range(0, 2**samples-1):
                rm = i;

                r_bin = dec2bin(i,samples);

                for j in range(1,samples-1):
                    r = bin2dec( np.roll(r_bin,-1) )
                    if r < rm:
                        rm = r;
                if tmpMap(rm+1) < 0:
                    tmpMap(rm+1) = newMax;
                    newMax = newMax + 1;
                table(i+1) = tmpMap(rm+1);

        if mappingtype == 'riu2': #Uniform & Rotation invariant
            newMax = samples + 2;
            for i in range(0, 2**samples - 1):

                i_bin =  dec2bin(i,samples);
                j_bin = np.roll(i_bin,-1);
                numt = sum(i_bin ^ j_bin);

                if numt <= 2:
                    table(i+1) = sum (bitget(i, 1:samples) ); #### TODO
                else:
                    table(i+1) = samples + 1;


########################################################## LBP.m

def lbp(image,radius=None, neighbors=8, mapping=0, mode='h'):
    spoints=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];

    if radius != None:
        radius=varargin{2};
        neighbors=varargin{3};

        spoints=np.zeros(neighbors, 2);

        # Angle step.
        a = 2*np.pi / neighbors;

        for i in range(1:neighbors):
            spoints(i,1) = -radius*sin((i-1)*a);
            spoints(i,2) = radius*cos((i-1)*a);

    ysize, xsize, depth = image.shape
