#! /usr/bin/python

import sys
import struct
from PyQt4.QtCore import *
from PyQt4.QtGui import *

__author__="remis"
__date__ ="$17-Feb-2010 21:21:31$"

class Converter(QDialog):
    """
    Convert list of intervals into scn format file.
    """

    def __init__(self, parent=None):
        super(Converter, self).__init__(parent)

        self.iscanver = -103
        self.amplitudes = []
        self.intervals = []
        self.options = []

        self.to_filename = ''
        self.nint = None
        self.calfac2 = 1.0
        self.ioffset = 154
        self.ffilt = -1.0
	self.rms = 0.0
	self.treso = 0.0
        self.tresg = 0.0
        self.title = "Converted from text file containing event list: amplitudes, intervals."   #als."    #character*70
	self.expdate = '00-ooo-0000'    #character*11
	self.tapeID = 'Converted list          '      #character*24
	self.ipatch = 0              #integer32
	self.Emem = 0.0              #float
	self.cjump = 0               #integer32
	self.avamp = 1.0             #float

        grid = QGridLayout()
        grid.addWidget(QLabel("From file:"), 0, 0)
        grid.addWidget(QLabel("To file:"), 1, 0)
        self.txt_from_file = QLineEdit()
        grid.addWidget(self.txt_from_file, 0, 1)
        self.txt_to_file = QLineEdit()
        grid.addWidget(self.txt_to_file, 1, 1)
        but_from_file = QPushButton("Browse")
        grid.addWidget(but_from_file, 0, 2)
        but_to_file = QPushButton("Browse")
        grid.addWidget(but_to_file, 1, 2)
        but_convert = QPushButton("Convert")
        grid.addWidget(but_convert, 2, 1)
        but_quit = QPushButton("Quit")
        grid.addWidget(but_quit, 2, 2)
        self.setLayout(grid)

        self.connect(but_quit, SIGNAL("clicked()"), self, SLOT("close()"))
        self.connect(but_convert, SIGNAL("clicked()"), self.convert)
        self.connect(but_from_file, SIGNAL("clicked()"), self.browse_from_file)
        self.connect(but_to_file, SIGNAL("clicked()"), self.browse_to_file)

        self.setWindowTitle("Convert to SCN file...")

    def browse_from_file(self):
        ""
        filename = QFileDialog.getOpenFileName(self,
                "Open a text file...", ".",
                "TXT files (*.txt);;All files (*.*)")
        if filename:
            self.txt_from_file.setText(filename)
            self.load_from_DO_file(filename)
            #self.fill_series_list(self.data.series_names())

    def browse_to_file(self):
        ""
        self.to_filename = QFileDialog.getSaveFileName(self,
                "Save as SCN file...", ".scn",
                "SCN files (*.scn)")
        if self.to_filename:
            #self.data.load_from_txt_file(filename)
            #self.fill_series_list(self.data.series_names())
            self.txt_to_file.setText(self.to_filename)

    def convert(self):
        ""
        self.save_binary()

    def save_binary(self):
        ""

        fout = open(self.to_filename, 'wb')
        fout.write(struct.pack('i', self.iscanver))
        fout.write(struct.pack('i', self.ioffset))
        fout.write(struct.pack('i', self.nint))

        fout.write(self.title)
        fout.write(self.expdate)
        fout.write(self.tapeID)

        fout.write(struct.pack('i', self.ipatch))
        fout.write(struct.pack('f', self.Emem))
        fout.write(struct.pack('i', 0))
        fout.write(struct.pack('f', self.avamp))
        fout.write(struct.pack('f', self.rms))
        fout.write(struct.pack('f', self.ffilt))
        fout.write(struct.pack('f', self.calfac2))
        fout.write(struct.pack('f', self.treso))
        fout.write(struct.pack('f', self.tresg))

#        for i in range(0, self.nint):
#            fout.write(struct.pack('f', self.intervals[i]))
#        for i in range(0, self.nint):
#            fout.write(struct.pack('h', self.amplitudes[i]))
#        for i in range(0, self.nint):
#            fout.write(struct.pack('b', 0))

        for i in range(0, self.nint):
            fout.write(struct.pack('f', self.intervals[i]))
        for i in range(0, self.nint):
            fout.write(struct.pack('h', self.amplitudes[i]))
        for i in range(0, self.nint):
            fout.write(struct.pack('b', 0))

        #currentff.byteswap()
        #currentff.tofile(fout)

        fout.close()
        print 'Finished writing scn file.'


    def save_qt_binary(self):
        error = None
        fh = None
        try:
            fh = QFile(self.to_filename)
            if not fh.open(QIODevice.WriteOnly):
                raise IOError, unicode(fh.errorString())
            stream = QDataStream(fh)
            stream.setByteOrder(QDataStream.LittleEndian)
            stream.writeInt32(self.iscanver)
            stream.writeInt32(self.ioffset)
            #stream.setVersion(QDataStream.Qt_4_2)
            stream.writeInt32(self.nint)
            for i in range(0, 105):
                stream.writeInt8('0')
            #title = str(self.title)
            #stream << title
            #stream.writeString(str(self.title))
            #stream.writeString(str(self.expdate))
            #stream.writeString(str(self.tapeID))
            stream.writeInt32(self.ipatch)
            stream.writeFloat(self.Emem)
            stream.writeInt32(0)
            #stream.writeInt32(self.cjump)
            stream.writeFloat(self.avamp)
            stream.writeFloat(self.rms)
            stream.writeFloat(self.ffilt)
            stream.writeFloat(self.calfac2)
            stream.writeFloat(self.treso)
            stream.writeFloat(self.tresg)
            for i in range(0, self.nint):
                stream.writeFloat(self.intervals[i])
            for i in range(0, self.nint):
                stream.writeInt16(self.amplitudes[i])
            for i in range(0, self.nint):
                stream.writeInt8(self.options[i])

        except (IOError, OSError), e:
            error = "Failed to save: {0}".format(e)
        finally:
            if fh is not None:
                fh.close()
                print 'saved'
            if error is not None:
                return False, error

        self.read_scn_file()


    def read_scn_file(self):
        ""

        error = None
        fh = None
        try:
            fh = QFile(self.to_filename)
            if not fh.open(QIODevice.ReadOnly):
                raise IOError, unicode(fh.errorString())
            stream = QDataStream(fh)
            stream.setByteOrder(QDataStream.LittleEndian)
            version = stream.readInt32()
            if version != self.iscanver:
                raise IOError, "version is wrong"
            else:
                print 'version =', version
            offset = stream.readInt32()
            if offset != self.ioffset:
                raise IOError, "offset is wrong"
            else:
                print 'offset =', offset
            nint = stream.readInt32()
            print 'number of intervals =', nint
            #stream.setVersion(QDataStream.Qt_4_2)

            #print stream.readInt32()
            print stream.readString()
            print stream.readString()
            print stream.readString()
            stream.readInt32()
            stream.readFloat()
            stream.readInt32()
            stream.readFloat()
            stream.readFloat()
            stream.readFloat()
            stream.readFloat()
            stream.readFloat()
            stream.readFloat()

            for i in range(self.nint):
                print (stream.readFloat())
                print (stream.readInt16())
                print (stream.readInt8())


        except (IOError, OSError), e:
            error = "Failed to load: {0}".format(e)
        finally:
            if fh is not None:
                fh.close()
                print 'scn file read well'
            if error is not None:
                return False, error


    def load_from_txt_file(self, filename=None):
        'Read data from a text file in which amplitudes and corresponding \
        interval length are in columns and are Tab delimited.'

        f = open(filename, 'r')
        lines =f.readlines()
        f.close()
        self.amplitudes = []
        self.intervals = []
        self.options = []
        for line in lines:
            line = line.strip("\n")    #remove newline
            values=line.split('\t')    #divide lines into values at tabs
            self.amplitudes.append(int(values[0]))
            self.intervals.append(float(values[1]))
            self.options.append(0)

        self.nint = len(self.intervals)
        print "number of imported intervals =", self.nint

#        for i in range(0, 19):
#            print self.amplitudes[i], self.intervals[i], self.options[i]

    def load_from_DO_file(self, filename=None):
        'Read data from a text file in which amplitudes and corresponding \
        interval length are in columns and are Tab delimited.'

        f = open(filename, 'r')
#        lin = f.readline()
        lines = f.readlines()
        f.close()
        self.amplitudes = []
        self.intervals = []
        self.options = []
        count = 0
        for line in lines:
            li = line.strip("\n")    #remove newline
            value = float(li.strip("\t"))    #divide lines into values at tabs
            print 'value=', value
            count += 1

            self.amplitudes.append(1)
            self.intervals.append(100)
            self.options.append(struct.pack("b",0)) ## int8_t


            self.amplitudes.append(0)
            self.intervals.append(value)
            self.options.append(struct.pack("=b",0))

        self.amplitudes.append(1)
        self.intervals.append(100)
        self.options.append(struct.pack("=b",0))

        self.nint = len(self.intervals)
        print "number of imported intervals =", self.nint
        print "number of useful intervals =", count


if __name__ == "__main__":
    print "Hello World"
    app = QApplication(sys.argv)
    converter = Converter()
    converter.show()
    app.exec_()
