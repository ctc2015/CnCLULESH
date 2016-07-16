
//********************************************************************************
// Copyright (c) 2010-2014 Intel Corporation. All rights reserved.              **
//                                                                              **
// Redistribution and use in source and binary forms, with or without           **
// modification, are permitted provided that the following conditions are met:  **
//   * Redistributions of source code must retain the above copyright notice,   **
//     this list of conditions and the following disclaimer.                    **
//   * Redistributions in binary form must reproduce the above copyright        **
//     notice, this list of conditions and the following disclaimer in the      **
//     documentation and/or other materials provided with the distribution.     **
//   * Neither the name of Intel Corporation nor the names of its contributors  **
//     may be used to endorse or promote products derived from this software    **
//     without specific prior written permission.                               **
//                                                                              **
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"  **
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE    **
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE   **
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE     **
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR          **
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF         **
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS     **
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      **
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)      **
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF       **
// THE POSSIBILITY OF SUCH DAMAGE.                                              **
//********************************************************************************
//

#ifndef CNC_COMMON_H
#define CNC_COMMON_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cnc/cnc.h>
#include <cassert>
#include <cnc/debug.h>
#include <memory>

#define TOI3D( _i, _j, _k, _szj, _szk ) ((_i)*(_szj)*(_szk)+(_j)*(_szk)+(_k))
#define TOI2D( _i, _j, _s ) ((_i)*(_s)+(_j))


typedef std::pair<int, int> Pair;
// don't use a vector: tags are copied and vector-copy is expensive
class Triple
{
public:
    Triple() { m_arr[0] = m_arr[1] = m_arr[2] = 0; }
    Triple( int a, int b, int c ) { m_arr[0] = a; m_arr[1] = b; m_arr[2] = c; }
    inline int operator[]( int i )  const { return m_arr[i]; }
private:
    int m_arr[3];
};


template <>
class cnc_tag_hash_compare< Triple >
{
public:
    size_t hash( const Triple & t ) const
    {
        return (t[0] << 11) + (t[1]) + ( t[2] << 17);      
    }
    bool equal( const Triple & a, const Triple & b) const { return a[0] == b[0] && a[1] == b[1] && a[2] == b[2]; }
};

inline std::ostream & cnc_format( std::ostream& os, const Triple & t )
{
    os << "[ " << t[0] << ", " << t[1] << ", " << t[2] << " ]";
    return os;
}

class Quad
{
public:
    Quad() { m_arr[0] = m_arr[1] = m_arr[2] = m_arr[3] = 0; }
    Quad( int a, int b, int c, int d) { m_arr[0] = a; m_arr[1] = b; m_arr[2] = c; m_arr[3] = d;}
    inline int operator[]( int i )  const { return m_arr[i]; }
private:
    int m_arr[4];
};


template <>
class cnc_tag_hash_compare< Quad >
{
public:
    size_t hash( const Quad & t ) const
    {
        return (t[0] << 16) + (t[1] << 8) + ( t[2] << 24) + t[3];      
    }
    bool equal( const Quad & a, const Quad & b) const { return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3]; }
};

inline std::ostream & cnc_format( std::ostream& os, const Quad & q )
{
    os << "[ " << q[0] << ", " << q[1] << ", " << q[2] << ", " << q[3] << " ]";
    return os;
}

template< typename T >
class Tile3d
{
public:
    Tile3d( int sz_x = 0, int sz_y = 0, int sz_z = 0) 
    : m_sz_x( sz_x ),
      m_sz_y( sz_y ),
      m_sz_z( sz_z ),
      m_array( 0 )
    { if( sz_x && sz_y && sz_z) m_array.resize( sz_x*sz_y*sz_z ); }

    ~Tile3d() {}
    inline T operator()( int i, int j, int k ) const { 
        return m_array[TOI3D(i,j,k,m_sz_y,m_sz_z)]; 
    }
    inline T & operator()( int i, int j, int k ) { 
        return m_array[TOI3D(i,j,k,m_sz_y,m_sz_z)]; 
    }

private:
    Tile3d( const Tile3d< T > & ) { assert( 0 ); }
    Tile3d & operator=( const Tile3d< T > & ) { assert( 0 ); return *this; }
    
    int   m_sz_x;
    int   m_sz_y;
    int   m_sz_z;
    std::vector< T > m_array;
    //    bool  m_full;
public:
};

#endif /* CNC_COMMON_H */
