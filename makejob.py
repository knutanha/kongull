from string import Template
textfile = open('templatejob', 'r')
template = textfile.read()

#print template

# data = [MPI, OMP, SIZE, MIN, SEK, LNODES, PPN, part]

data = [
[1, 1, 2048,0,20,1,1,'B'],
[3, 12, 2048,0,20,3,12,'B'],
[6, 6, 2048,0,20,3,12,'B'],
[36, 1, 2048,0,20,3,12,'B'],

[3, 12, 16384,7,0,3,12,'C'],
[6, 6, 16384,7,0,3,12,'C'],
[12, 3, 16384,7,0,3,12,'C'],
[36, 1, 16384,7,0,3,12,'C'],

[1, 1, 16384,10,0,1,1,'Dsp'],
[3, 1, 16384,7,0,3,1,'Dsp'],
[6, 1, 16384,7,0,3,2,'Dsp'],
[6, 2, 16384,7,0,3,4,'Dsp'],
[6, 4, 16384,7,0,3,8,'Dsp'],
#[6, 6, 16384,0,20,3,12,'D']

[6, 6, 128,0,15,3,12,'Dsn'],
[6, 6, 256,0,15,3,12,'Dsn'],
[6, 6, 512,0,15,3,12,'Dsn'],
[6, 6, 1024,0,20,3,12,'Dsn'],
[6, 6, 2048,0,30,3,12,'Dsn'],
[6, 6, 4096,2,0,3,12,'Dsn'],
[6, 6, 8192,4,0,3,12,'Dsn']
#[6, 6, 16384,0,20,3,12,'D']
]

for i in xrange(0,len(data)):
    name = "%s_%i_%i_%i" % (data[i][7], data[i][0], data[i][1], data[i][2])
    mpipn = "%i" % (data[i][0]/data[i][5])
    omp = "%i" % (data[i][1])
    size = "%i" % (data[i][2])
    minutes = "%i" % (data[i][3])
    secs = "%i" % (data[i][4])
    lnodes = "%i" % (data[i][5])
    ppn = "%i" % (data[i][6])
    test = Template(template)
    a = test.safe_substitute(NAME=name,MIN=minutes, SEK=secs, LNODES=lnodes, PPN=ppn, OMP=omp,MPIPN=mpipn,SIZE=size)
    name = name + '.sh'
    newfile = open(name,'w')
    newfile.write(a)
    newfile.close()

