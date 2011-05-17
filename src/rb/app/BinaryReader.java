//    rbAPPmit: An Android front-end for the Certified Reduced Basis Method
//    Copyright (C) 2010 David J. Knezevic and Phuong Huynh
//
//    This file is part of rbAPPmit
//
//    rbAPPmit is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rbAPPmit is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rbAPPmit.  If not, see <http://www.gnu.org/licenses/>. 

package rb.app;

import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;

/* source: http://www.captain.at/howto-java-convert-binary-data.php */

public class BinaryReader extends DataInputStream{
	
	public BinaryReader(InputStream in) {
		super(in);
		// TODO Auto-generated constructor stub
	}
	
	// Read a float number
	public float ReadFloat(){
		int i = 0;
		byte[] tmp = new byte[4];
		try {
			for (i = 0; i < 4; i++)
				tmp[i] = readByte();
		} catch (IOException e) {
			throw new RuntimeException("Error reading data", e);
		}		
		int accum = 0;
		i = 0;
		for (int shiftBy = 0; shiftBy < 32; shiftBy += 8 ) {
			accum |= ( (long)( tmp[i] & 0xff ) ) << shiftBy;
			i++;
		}
		return Float.intBitsToFloat(accum);
	}
	
	// Read a float array
	public float[] ReadFloat(int _size){
		float[] ofloat = new float[_size];
		int i = 0;
		byte[] tmp = new byte[4*_size];
		try {
			for (i = 0; i < 4*_size; i++)
				tmp[i] = readByte();
		} catch (IOException e) {
			throw new RuntimeException("Error reading data", e);
		}
		int accum;
		for (int count = 0; count < _size; count++){
			accum = 0;
			i = 0;
			for (int shiftBy = 0; shiftBy < 32; shiftBy += 8 ) {
				accum |= ( (long)( tmp[i+count*4] & 0xff ) ) << shiftBy;
				i++;
			}
			ofloat[count] = Float.intBitsToFloat(accum);	
		}
		return ofloat;		
	}
	
	// Read a double number
	public double ReadDouble(){
		int i = 0;
		byte[] tmp = new byte[8];
		try {
			for (i = 0; i < 8; i++)
				tmp[i] = readByte();
		} catch (IOException e) {
			throw new RuntimeException("Error reading data", e);
		}
		long accum = 0;
		i = 0;
		for ( int shiftBy = 0; shiftBy < 64; shiftBy += 8 ) {
			accum |= ( (long)( tmp[i] & 0xff ) ) << shiftBy;
			i++;
		}
		return Double.longBitsToDouble(accum);
	}
	
	// Read an integer number
	public int ReadInt(){
		int i = 0;
		byte[] tmp = new byte[2];
		try {
			for (i = 0; i < 2; i++)
				tmp[i] = readByte();
		} catch (IOException e) {
			throw new RuntimeException("Error reading data", e);
		}
		int low = tmp[0] & 0xff;
		int high = tmp[1] & 0xff;
		return (int)( high << 8 | low );
	}
	
	// Read a long number
	public long ReadLong(){
		int i = 0;
		byte[] tmp = new byte[4];
		try {
			for (i = 0; i < 4; i++)
				tmp[i] = readByte();
		} catch (IOException e) {
			throw new RuntimeException("Error reading data", e);
		}
		long accum = 0;
		i = 0;
		for ( int shiftBy = 0; shiftBy < 32; shiftBy += 8 ) {
			accum |= ( (long)( tmp[i] & 0xff ) ) << shiftBy;
			i++;
		}
		return accum;
	}
}

