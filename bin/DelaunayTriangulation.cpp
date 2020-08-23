#include "DelaunayTriangulation.h"

DelaunayTriangulation::DelaunayTriangulation(double u0, double u1, double v0, double v1, NURBSSurface surface, double eps, double metr) {
	this->id = 1;
	this->surface = surface;
	this->eps = eps;
	this->metr = metr;
	std::array<double, 2> a = { u0 - 1, v0 - 1 };
	std::array<double, 2> b = { u1 + 1, v0 - 1 };
	std::array<double, 2> c = { u1 + 1, v1 + 1 };
	std::array<double, 2> d = { u0 - 1, v1 + 1 };
	this->coords.push_back(a);
	this->coords.push_back(b);
	this->coords.push_back(c);
	this->coords.push_back(d);
	std::array<int, 3> T1 = { 0, 1 ,3 };
	std::array<int, 3> T2 = { 2, 3, 1 };
	this->triangles.emplace(0, T1);
	this->triangles.emplace(1, T2);
	std::array<int, 3> bound_T1 = { 1, -1, -1 };
	std::array<int, 3> bound_T2 = { 0, -1, -1 };
	this->neighbor_triangles.emplace(0, bound_T1);
	this->neighbor_triangles.emplace(1, bound_T2);
	this->circles.emplace(0, circumcenter(T1));
	this->circles.emplace(1, circumcenter(T2));
}

std::array<double, 3> DelaunayTriangulation::circumcenter(std::array<int, 3> tri) {
	std::array<double, 2> a = this->coords.at(tri[0]);
	std::array<double, 2> b = this->coords.at(tri[1]);
	std::array<double, 2> c = this->coords.at(tri[2]);
	double d = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
	std::array<double, 3> center = { ((a[0] * a[0] + a[1] * a[1]) * (b[1] - c[1]) + (b[0] * b[0] + b[1] * b[1]) * (c[1] - a[1]) + (c[0] * c[0] + c[1] * c[1]) * (a[1] - b[1])) / d,
						((a[0] * a[0] + a[1] * a[1]) * (c[0] - b[0]) + (b[0] * b[0] + b[1] * b[1]) * (a[0] - c[0]) + (c[0] * c[0] + c[1] * c[1]) * (b[0] - a[0])) / d, 0 };
	center[2] = (a[0] - center[0]) * (a[0] - center[0]) + (a[1] - center[1]) * (a[1] - center[1]);
	return center;
}

bool DelaunayTriangulation::inCircle(int tri_id, std::array<double, 2> p) {
	double eps = 0.0000001;
	std::array<double, 3> center = this->circles.at(tri_id);
	return ((p[0] - center[0]) * (p[0] - center[0]) + (p[1] - center[1]) * (p[1] - center[1])) <= center[2] + eps;
}

void DelaunayTriangulation::addPoint(std::array<double, 2> p, bool euclid_flag) {
	int idx = this->coords.size();
	this->coords.push_back(p);
	std::vector<int> bad_triangles;
	int T_id;
	if (idx == 234) {
		for (auto i : this->coords) { std::cout << "(" << i[0] << ", " << i[1] << ")," << std::endl; }
		std::cout << std::endl << std::endl << std::endl;
		for (auto i : this->mesh) {
			std::cout << "(" << i.second[0] << ", " << i.second[1] << ", " << i.second[2] << ")," << std::endl;
		}
		std::cout << std::endl << std::endl << std::endl;
		for (auto i : this->neighbor_triangles) {
			std::cout << "(" << this->triangles[i.first][0] << ", " << this->triangles[i.first][1] << ", " << this->triangles[i.first][2] << ")," << std::endl;
		}
	}
	//std::cout << std::endl << std::endl << std::endl;
	//for (auto i : this->coords) { std::cout << "(" << i[0] << ", " << i[1] << ")," << std::endl; }
	//std::cout << std::endl << std::endl << std::endl;
	for (auto tri : this->neighbor_triangles) {
		if (this->initTriangle(tri.first, p)) {
			T_id = tri.first;
			break;
		}
	}
	bad_triangles.push_back(T_id);
	std::vector<std::array<int, 3>> boundary;
	int edge = 0;
	int tri_op_id;
	int prev = 0;
	int next = 0;
	int pos = 0;
	//for (auto i : bad_triangles) {
	//	std::cout << "(" << this->triangles[i][0] << ", " << this->triangles[i][1] << ", " << this->triangles[i][2] << ")," << std::endl;
	//}
	//std::cout << std::endl << std::endl << std::endl;
	if (euclid_flag) {
		while (true) {
			tri_op_id = this->neighbor_triangles[T_id][edge];
			if (tri_op_id==-1 || !this->inCircle(tri_op_id, p)) {
				if (edge == 0) { prev = 2; next = 1; }
				if (edge == 1) { prev = 0; next = 2; }
				if (edge == 2) { prev = 1; next = 0; }
				std::array<int, 3> item = { this->triangles[T_id][next], this->triangles[T_id][prev], tri_op_id };
				boundary.push_back(item);
				edge = (edge + 1) % 3;
				if (boundary[0][0] == boundary[boundary.size() - 1][1]) { break; }
			}
			else {
				if (!this->checkVector(tri_op_id, bad_triangles)) { bad_triangles.push_back(tri_op_id); }
				for (int i = 0; i < 3; i++) {
					if (T_id == this->neighbor_triangles[tri_op_id][i]) {
						pos = i;
						break;
					}
				}
				edge = (pos + 1) % 3;
				T_id = tri_op_id;
			}
		}
	}
	else {
		while (true) {
			tri_op_id = this->neighbor_triangles[T_id][edge];
			if (tri_op_id == -1 || !this->inCircumdisk(tri_op_id, p)) {
				if (edge == 0) { prev = 2; next = 1; }
				if (edge == 1) { prev = 0; next = 2; }
				if (edge == 2) { prev = 1; next = 0; }
				std::array<int, 3> item = { this->triangles[T_id][next], this->triangles[T_id][prev], tri_op_id };
				boundary.push_back(item);
				edge = (edge + 1) % 3;
				if (boundary[0][0] == boundary[boundary.size() - 1][1]) { break; }
			}
			else {
				if (!this->checkVector(tri_op_id, bad_triangles)) { bad_triangles.push_back(tri_op_id); }
				for (int i = 0; i < 3; i++) {
					if (T_id == this->neighbor_triangles[tri_op_id][i]) {
						pos = i;
						break;
					}
				}
				edge = (pos + 1) % 3;
				T_id = tri_op_id;
			}
		}
	}
	int count = 0;
	int step = -1;
	std::vector<int> new_triangles;
	for (auto item : boundary) {
		std::array<int, 3> T = { idx, item[0], item[1] };
		this->id += 1;
		if (this->id == 1906) {
			std::cout << std::endl;
		}
		this->triangles.emplace(this->id, T);
		if(euclid_flag)
			this->circles.emplace(this->id, this->circumcenter(T));
		else
			this->circles.emplace(this->id, this->circumdiskCenter(T));
		std::array<int,3> n_T = { item[2], -1, -1 };
		this->neighbor_triangles.emplace(this->id, n_T);
		if (item[2] != -1) {
			for (auto i : this->neighbor_triangles[item[2]]) {
				step++;
				count = 0;
				if (i != -1) {
					for (int j = 0; j < 3; j++) {
						if ((this->triangles[i][j] == item[0]) || (this->triangles[i][j] == item[1]))
							count++;
					}
				}
				if (count == 2) {
					this->neighbor_triangles[item[2]][step] = this->id;
					step = -1;
					break;
				}
			}
		}
		new_triangles.push_back(this->id);
	}
	step = -1;
	for (auto T : bad_triangles) {
		this->delete_id.push_back(T);
		this->neighbor_triangles.erase(T);
		this->circles.erase(T);
	}
	for (auto T : new_triangles) {
		step++;
		prev = step - 1;
		next = (step + 1) % new_triangles.size();
		if (prev == -1) { prev = new_triangles.size() - 1; }
		this->neighbor_triangles[T][1] = new_triangles[next];
		this->neighbor_triangles[T][2] = new_triangles[prev];
	}
}

std::map<int, std::array<int, 3>> DelaunayTriangulation::exportTriangles() {
	return this->triangles;
}

void DelaunayTriangulation::boundaryRecovery(std::vector<std::array<double, 2>> front, bool frame_flag) {
	int id;
	std::vector<int> vec;
	for (int i : this->delete_id)
		this->triangles.erase(i);
	this->delete_id.clear();
	if (frame_flag) {
		for (auto i : this->triangles) {
			for (auto j : i.second) {
				if ((j == 0) || (j == 1) || (j == 2) || (j == 3)) {
					id = i.first;
					vec.push_back(id);
					break;
				}
			}
			for (auto k : this->neighbor_triangles) {
				for (int l = 0; l < 3; l++) {
					if (k.second[l] == id)
						this->neighbor_triangles[k.first][l] = -1;
				}
			}
		}
		for (int i : vec) {
			this->triangles.erase(i);
			this->neighbor_triangles.erase(i);
			this->circles.erase(i);
		}
	}
	else {
		std::array<double, 2> A;
		std::array<double, 2> B;
		int a_id;
		int b_id;
		int c_id;
		int count;
		for (int i=0; i < front.size() - 1; i++) {
			A = front[i];
			B = front[i + 1];
			count = 0;
			for (auto item : front) {
				if (item != A && item != B && CurveLength(A, item) < (this->metr + this->eps) && CurveLength(B, item) < (this->metr + this->eps)) {
					for (int i = 0; i < this->coords.size(); i++) {
						if (this->coords.at(i) == A)
							a_id = i;
						if (this->coords.at(i) == B)
							b_id = i;
						if (this->coords.at(i) == item)
							c_id = i;
					}
					std::array<int, 3> new_T = { a_id,b_id,c_id };
					this->mesh.emplace(this->mesh.size(), new_T);
					count++;
					if (count > 1) { break; }
				}
			}
		}
	}
}

bool DelaunayTriangulation::checkVector(int i, std::vector<int> vec) {
	for (auto j : vec) {
		if (j == i)
			return true;
	}
	return false;
}

double DelaunayTriangulation::ScalarMulV3(double* v1, double* v2) {
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

double DelaunayTriangulation::CurveLength(std::array<double,2> A, std::array<double, 2> B) {
	double l1[3];
	double l2[3];
	double r = 0;
	double dt = 0.1;
	double AB[2] = { B[0] - A[0], B[1] - A[1] };
	for (double t = 0; t < 1 - this->eps; t += dt) {
		this->surface.getDer1Point(A[0] + t * AB[0], A[1] + t * AB[1], l1);
		this->surface.getDer2Point(A[0] + t * AB[0], A[1] + t * AB[1], l2);
		double metric_A[2][2] = { {ScalarMulV3(l1,l1), ScalarMulV3(l1, l2)}, {ScalarMulV3(l2, l1), ScalarMulV3(l2, l2) } };
		r += pow((AB[0] * metric_A[0][0] + AB[1] * metric_A[1][0]) * AB[0] + (AB[0] * metric_A[0][1] + AB[1] * metric_A[1][1]) * AB[1], 0.5) * dt;
	}
	return r;
}

std::array<double, 3> DelaunayTriangulation::circumdiskCenter(std::array<int, 3> tri) {
	eps = 0.1;
	std::array<double, 2> a = this->coords.at(tri[0]);
	std::array<double, 2> b = this->coords.at(tri[1]);
	std::array<double, 2> c = this->coords.at(tri[2]);
	double d = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
	std::array<double, 2> center = { ((a[0] * a[0] + a[1] * a[1]) * (b[1] - c[1]) + (b[0] * b[0] + b[1] * b[1]) * (c[1] - a[1]) + (c[0] * c[0] + c[1] * c[1]) * (a[1] - b[1])) / d,
						((a[0] * a[0] + a[1] * a[1]) * (c[0] - b[0]) + (b[0] * b[0] + b[1] * b[1]) * (a[0] - c[0]) + (c[0] * c[0] + c[1] * c[1]) * (b[0] - a[0])) / d };
	
	while (1) {
		double da = this->CurveLength(center, a);
		double db = this->CurveLength(center, b);
		double dc = this->CurveLength(center, c);
		if (da > (db - eps) && da < (db + eps) && da >(dc - eps) && da < (dc + eps) && dc >(db - eps) && dc < (db + eps))
			break;
		double d = (da + db + dc) / 3;
		std::array<double, 2> ao = { d * (center[0] - a[0]) / da,d * (center[1] - a[1]) / da };
		std::array<double, 2> bo = { d * (center[0] - b[0]) / db,d * (center[1] - b[1]) / db };
		std::array<double, 2> co = { d * (center[0] - c[0]) / dc,d * (center[1] - c[1]) / dc };
		std::array<double, 2> oa = { a[0] + ao[0],a[1] + ao[1] };
		std::array<double, 2> ob = { b[0] + bo[0],b[1] + bo[1] };
		std::array<double, 2> oc = { c[0] + co[0],c[1] + co[1] };
		center[0] = (oa[0] + ob[0] + oc[0]) / 3;
		center[1] = (oa[1] + ob[1] + oc[1]) / 3;
	}

	std::array<double, 3> answer = { center[0], center[1], 0.0 };
	return answer;
}

bool DelaunayTriangulation::initTriangle(int tri_id, std::array<double, 2> p) {
	std::array<int, 3> tri = this->triangles.at(tri_id);
	std::array<double, 2> a = this->coords.at(tri[0]);
	std::array<double, 2> b = this->coords.at(tri[1]);
	std::array<double, 2> c = this->coords.at(tri[2]);
	double p1 = (a[0] - p[0]) * (b[1] - a[1]) - (b[0] - a[0]) * (a[1] - p[1]);
	double p2 = (b[0] - p[0]) * (c[1] - b[1]) - (c[0] - b[0]) * (b[1] - p[1]);
	double p3 = (c[0] - p[0]) * (a[1] - c[1]) - (a[0] - c[0]) * (c[1] - p[1]);
	if ((p1 >= 0 && p2 >= 0 && p3 >= 0) || (p1 <= 0 && p2 <= 0 && p3 <= 0)) { return true; }
	return false;
}

bool DelaunayTriangulation::inCircumdisk(int tri_id, std::array<double, 2> p) {
	double eps = 0.001;
	std::array<double, 3> center = this->circles.at(tri_id);
	std::array<double, 2> o = { center[0], center[1] };
	std::array<int, 3> tri = this->triangles.at(tri_id);
	std::array<double, 2> a = this->coords.at(tri[0]);
	std::array<double, 2> b = this->coords.at(tri[1]);
	std::array<double, 2> c = this->coords.at(tri[2]);

	double oa = this->CurveLength(o, a);
	double ob = this->CurveLength(o, b);
	double oc = this->CurveLength(o, c);
	double max = oa;
	if (ob > max) {
		max = ob;
	}
	if (oc > max) {
		max = oc;
	}

	return this->CurveLength(o, p) <= max + eps;
}

void DelaunayTriangulation::recalculateCircumcentre() {
	std::vector<int> indices;
	for (auto item : this->circles)
		indices.push_back(item.first);
	for (auto index : indices) {
		this->circles.erase(index);
		this->circles.emplace(index, circumdiskCenter(this->triangles[index]));
	}
}

void DelaunayTriangulation::deleteBaseTriangle(std::array<double, 2> a, std::array<double, 2> b, std::array<double, 2> c) {
	int a_id;
	int b_id;
	int c_id;
	for (int i = 0; i < this->coords.size(); i++) {
		if (this->coords.at(i) == a)
			a_id = i;
		if (this->coords.at(i) == b)
			b_id = i;
		if (this->coords.at(i) == c)
			c_id = i;
	}
	int t;
	std::array<int, 3> tri = { a_id, b_id, c_id };
	std::array<int, 3> tri1 = { c_id, a_id, b_id };
	std::array<int, 3> tri2 = { b_id, c_id, a_id };
	for (auto item : this->triangles) {
		if (item.second == tri || item.second == tri1 || item.second == tri2) {
			t = item.first;
			break;
		}
	}
	this->mesh.emplace(this->mesh.size(), this->triangles[t]);
	this->delete_id.push_back(t);
	this->neighbor_triangles.erase(t);
	this->circles.erase(t);
	for (auto k : this->neighbor_triangles) {
		for (int l = 0; l < 3; l++) {
			if (k.second[l] == t)
				this->neighbor_triangles[k.first][l] = -1;
		}
	}
}

std::vector<std::array<double, 2>> DelaunayTriangulation::getCoords() {
	return this->coords;
}

std::map<int, std::array<int, 3>> DelaunayTriangulation::getMesh() {
	return this->mesh;
}

std::array<double, 2> DelaunayTriangulation::selectFromDT(std::array<double, 2> a, std::array<double, 2> b) {
	int a_id;
	int b_id;
	int c_id;
	int tri_id;
	for (int i = 0; i < this->coords.size(); i++) {
		if (this->coords.at(i) == a)
			a_id = i;
		if (this->coords.at(i) == b)
			b_id = i;
	}
	for (auto item : this->neighbor_triangles) {
		if (this->triangles[item.first][0] == a_id && this->triangles[item.first][1] == b_id) {
			c_id = this->triangles[item.first][2];
			tri_id = item.first;
			break;
		}
		if (this->triangles[item.first][1] == a_id && this->triangles[item.first][2] == b_id) {
			c_id = this->triangles[item.first][0];
			tri_id = item.first;
			break;
		}
		if (this->triangles[item.first][2] == a_id && this->triangles[item.first][0] == b_id) {
			c_id = this->triangles[item.first][1];
			tri_id = item.first;
			break;
		}
	}
	this->mesh.emplace(this->mesh.size(), this->triangles[tri_id]);
	this->delete_id.push_back(tri_id);
	this->neighbor_triangles.erase(tri_id);
	this->circles.erase(tri_id);
	for (auto k : this->neighbor_triangles) {
		for (int l = 0; l < 3; l++) {
			if (k.second[l] == tri_id)
				this->neighbor_triangles[k.first][l] = -1;
		}
	}
	//for (auto i : this->mesh) {
	//	std::cout << "(" << i.second[0] << ", " << i.second[1] << ", " << i.second[2] << ")," << std::endl;
	//}
	//std::cout << std::endl << std::endl << std::endl;
	//for (auto i : this->neighbor_triangles) {
	//	std::cout << "(" << this->triangles[i.first][0] << ", " << this->triangles[i.first][1] << ", " << this->triangles[i.first][2] << ")," << std::endl;
	//}
	//std::cout << std::endl << std::endl << std::endl;
	//for (auto i : this->coords) { std::cout << "(" << i[0] << ", " << i[1] << ")," << std::endl; }
	return this->coords[c_id];
}

void DelaunayTriangulation::lastTriangle(std::array<double, 2> a, std::array<double, 2> b, std::array<double, 2> c) {
	int a_id;
	int b_id;
	int c_id;
	for (int i = 0; i < this->coords.size(); i++) {
		if (this->coords.at(i) == a)
			a_id = i;
		if (this->coords.at(i) == b)
			b_id = i;
		if (this->coords.at(i) == b)
			c_id = i;
	}
	std::array<int, 3> tri2 = { a_id, b_id, c_id };
	this->mesh.emplace(this->mesh.size(), tri2);
}

bool DelaunayTriangulation::actualPoints(std::array<double, 2> p) {
	this->act_points.clear();
	int ind;
	for (auto item : this->neighbor_triangles) {
		this->act_points.insert(this->triangles[item.first][0]);
		this->act_points.insert(this->triangles[item.first][1]);
		this->act_points.insert(this->triangles[item.first][2]);
	}
	for (int i = 0; i < this->coords.size(); i++) {
		if (this->coords[i] == p) {
			ind = i;
			break;
		}
	}
	if (this->act_points.find(ind) == this->act_points.end()) {
		return true;
	}
	return false;
}