// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - David Andrs
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef _ADAPT_PROJ_H_
#define _ADAPT_PROJ_H_

/// @defgroup hp-adapt hp-Adaptivity


/// Abstract class for projecting reference solution to the coarse mesh
///
/// NOTE: hex-specific
///
/// @ingroup hp-adapt
class Projection {
public:
	Projection(Solution *afn, Element *e, Shapeset *ss);
	virtual ~Projection();

	virtual double get_error(int split, int son, order3_t order) = 0;

protected:
	order3_t order;
	Mesh *mesh;
	Solution *sln;
	Element *base_elem;
	Quad3D *quad;
	Shapeset *ss;
	ShapeFunction *fu;
	ShapeFunction *fv;

	struct ProjItem {
		scalar coef;				// coef of the projection
		int idx;					// index of the shape function
	};

	ProjItem *vertex_proj;			// vertex projection
	ProjItem *edge_proj[12];		// edge projection
	ProjItem *face_proj[6];			// face projection
	ProjItem *bubble_proj;			// bubble projection

	ProjItem **proj;				// projection
	int proj_fns;					// number of funcions

	Trf *get_trf(int trf);
	void free_proj();
	void calc_projection(int split, int son, order3_t order);
	virtual void calc_vertex_proj(int split, int son) = 0;
	virtual void calc_edge_proj(int edge, int split, int son, order3_t order) = 0;
	virtual void calc_face_proj(int face, int split, int son, order3_t order) = 0;
	virtual void calc_bubble_proj(int split, int son, order3_t order) = 0;

	// FIXME: Hex-specific
	static const int NUM_TRF = 27;		// number of all possible transformations
	// idx of sons to visit
	static int vtx_son[NUM_TRF][8];
	static int edge_son[NUM_TRF][Hex::NUM_EDGES][2];
	static int face_son[NUM_TRF][Hex::NUM_FACES][4];
	static int int_son[NUM_TRF][8];
	// number of sons to go over
	static int edge_ns[8][Hex::NUM_EDGES];
	static int face_ns[8][Hex::NUM_FACES];
	static int int_ns[8];
	// transformations to apply
	static int edge_trf[8][Hex::NUM_EDGES][2];
	static int face_trf[8][Hex::NUM_FACES][4];
	static int int_trf[8][8];
};

#endif
