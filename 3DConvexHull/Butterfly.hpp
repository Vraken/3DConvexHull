#pragma once
//#include "Graph.hpp"
#include <map>

class Butterfly
{
private:
	Graph* graph;
public:
	Butterfly()
	{
		graph = nullptr;
	}

	Butterfly(Graph* g)
	{
		graph = g;
	}

	void Subdivise()
	{
		std::map<Edge*, Summit*> edge_vertex;

		std::vector<Summit*>* new_summits = new std::vector<Summit*>();
		std::vector<Edge*>* new_edges = new std::vector<Edge*>();
		std::vector<Face*>* new_faces = new std::vector<Face*>();

		std::vector<Summit*>* old_summits = graph->getSummitList();
		std::vector<Edge*>* old_edges = graph->getEdgeList();
		std::vector<Face*>* old_faces = graph->getFaceList();

		/*std::vector<Summit*>* old_summits = new std::vector<Summit*>(*graph->getSummitList());
		for (int i = 0; i < graph->getSummitList()->size(); ++i)
		{
			old_summits->push_back(new Summit(*graph->getSummitList()->at(i)));
		}*/

		/*std::vector<Edge*>* old_edges = new std::vector<Edge*>(*graph->getEdgeList());
		for (int i = 0; i < graph->getEdgeList()->size(); ++i)
		{
			old_edges->push_back(new Edge(*graph->getEdgeList()->at(i)));
		}*/

		/*std::vector<Face*>* old_faces = new std::vector<Face*>(*graph->getFaceList());
		for (int i = 0; i < graph->getFaceList()->size(); ++i)
		{
			old_faces->push_back(new Face(*graph->getFaceList()->at(i)));
		}*/

		//Creation of edge vertexes
		for (int i = 0; i < old_edges->size(); ++i)
		{
			switch (GetButterflyScheme(old_edges->at(i)))
			{
			case 0:
				edge_vertex[old_edges->at(i)] = GetBoundaryBarycentre(old_edges->at(i));
				break;
			case 1:
				edge_vertex[old_edges->at(i)] = GetFullSchemeBarycentre(old_edges->at(i));
				break;
			case 2:
				edge_vertex[old_edges->at(i)] = GetSemiBoundaryBarycentre(old_edges->at(i));
				break;
			case 3:
				edge_vertex[old_edges->at(i)] = GetSemiBoundaryBarycentre(old_edges->at(i));
				break;
			case 4:
				edge_vertex[old_edges->at(i)] = GetPyramidalBarycentre(old_edges->at(i));
				break;
			case 5:
				edge_vertex[old_edges->at(i)] = GetUnbalancedBarycentre(old_edges->at(i));
				break;
			}
		}

		//Creation of new edges between previous points and edge vertexes
		Summit* cur_summit;
		for (int i = 0; i < old_summits->size(); ++i)
		{
			new_summits->push_back(new Summit(*old_summits->at(i)));
			cur_summit = new_summits->back();
			cur_summit->setEdgesConnected(new std::vector<Edge*>());
			cur_summit->setFacesConnected(new std::vector<Face*>());

			for (int j = 0; j < old_summits->at(i)->getEdgesConnected()->size(); ++j)
			{
				if (std::find(new_summits->begin(), new_summits->end(), edge_vertex[old_summits->at(i)->getEdgesConnected()->at(j)]) == new_summits->end())
				{
					new_summits->push_back(edge_vertex[old_summits->at(i)->getEdgesConnected()->at(j)]);
					new_summits->back()->setEdgesConnected(new std::vector<Edge*>());
					new_summits->back()->setFacesConnected(new std::vector<Face*>());
				}

				std::vector<Summit*> *summits = new std::vector<Summit*>(
				{ cur_summit , edge_vertex[old_summits->at(i)->getEdgesConnected()->at(j)] });

				new_edges->push_back(new Edge(summits));
				new_edges->back()->setSummitsConnected(summits);
				new_edges->back()->setFacesConnected(new std::vector<Face*>());

				cur_summit->getEdgesConnected()->push_back(new_edges->back());
				edge_vertex[old_summits->at(i)->getEdgesConnected()->at(j)]->getEdgesConnected()->push_back(new_edges->back());
			}
		}

		//Creation of new edges between edge vertexes
		std::map<Edge*, Summit*>::iterator it;
		for (it = edge_vertex.begin(); it != edge_vertex.end(); ++it)
		{
			for (int i = 0; i < it->first->getFacesConnected()->size(); ++i)
			{
				for (int j = 0; j < it->first->getFacesConnected()->at(i)->getEdgesConnected()->size(); ++j)
				{
					if (it->first->getFacesConnected()->at(i)->getEdgesConnected()->at(j) == it->first)
					{
						continue;
					}

					std::vector<Summit*> *summits = new std::vector<Summit*>(
					{ it->second , edge_vertex[it->first->getFacesConnected()->at(i)->getEdgesConnected()->at(j)] });

					Edge* new_edge = new Edge(summits);

					if (EdgeExists(new_edge, new_edges) == nullptr)
					{
						new_edges->push_back(new_edge);
						new_edge->setSummitsConnected(summits);
						new_edge->setFacesConnected(new std::vector<Face*>());
							
						it->second->getEdgesConnected()->push_back(new_edge);
						edge_vertex[it->first->getFacesConnected()->at(i)->getEdgesConnected()->at(j)]->getEdgesConnected()->push_back(new_edge);
					}
				}
			}
		}

		//Creation of new faces with new edges
		std::map<Edge*, Summit*>::iterator it2;
		for (it2 = edge_vertex.begin(); it2 != edge_vertex.end(); ++it2)
		{
			Summit* s1 = it2->second;

			for (int i = 0; i < s1->getEdgesConnected()->size(); ++i)
			{
				Edge* e1 = s1->getEdgesConnected()->at(i);
				Summit* s2 = GetOppositeSummit(s1, e1);

				for (int j = 0; j < s1->getEdgesConnected()->size(); ++j)
				{
					if (i == j)
					{
						continue;
					}

					Edge* e2 = s1->getEdgesConnected()->at(j);
					Summit* s3 = GetOppositeSummit(s1, e2);

					Edge* e3 = FindEdgeFromSummits(s2, s3, new_edges);

					if (e3 == nullptr)
					{
						continue;
					}

					std::vector<Summit*>* summits = new std::vector<Summit*>({ s1, s2, s3 });
					std::vector<Edge*>* edges = new std::vector<Edge*>({ e1, e2, e3 });

					Face* f = new Face(summits, edges);

					if (FindFaceFromSummits(s1, s2, s3, new_faces) != nullptr)
					{
						continue;
					}

					new_faces->push_back(f);

					s1->getFacesConnected()->push_back(f);
					s2->getFacesConnected()->push_back(f);
					s3->getFacesConnected()->push_back(f);

					e1->getFacesConnected()->push_back(f);
					e2->getFacesConnected()->push_back(f);
					e3->getFacesConnected()->push_back(f);
				}
			}
		}

		graph->setSummitList(new_summits);
		graph->setEdgeList(new_edges);
		graph->setFaceList(new_faces);
	}

	Summit* GetOppositeSummit(Edge* e, Face* f)
	{
		for (int i = 0; i < f->getSummitsConnected()->size(); ++i)
		{
			if (f->getSummitsConnected()->at(i) != e->getSummitsConnected()->at(0)
				&& f->getSummitsConnected()->at(i) != e->getSummitsConnected()->at(1))
			{
				return f->getSummitsConnected()->at(i);
			}
		}

		return nullptr;
	}

	Summit* GetOppositeSummit(Summit* s, Edge* e)
	{
		if (e->getSummitsConnected()->at(0) != s)
		{
			return e->getSummitsConnected()->at(0);
		}
		return e->getSummitsConnected()->at(1);
	}

	int GetButterflyScheme(Edge* e)
	{
		/*if (e->getFacesConnected()->size() == 1)
		{
			return 0;
		}*/
		//if
		//{
			int nb[2] = { 0, 0 };
			for (int i = 0; i < e->getFacesConnected()->size(); ++i)		//for each face
			{
				for (int j = 0; j < e->getFacesConnected()->at(i)->getEdgesConnected()->size(); ++j)		//for each edge of the face
				{
					if (e->getFacesConnected()->at(i)->getEdgesConnected()->at(j)->getFacesConnected()->size() == 2)
					{
						++nb[i];
					}
				}
			}

			if (nb[0] == nb[1])
			{
				switch (nb[0])
				{
				case 3:
					return 1;
				case 2:
					return 2;
				case 1:
					return 2;
				}
			}
			else
			{
				if (nb[0] == 0 || nb[1] == 0)
				{
					return 0;
				}
				else if ((nb[0] == 1 && nb[1] == 2)
					|| (nb[0] == 2 && nb[1] == 1))
				{
					return 3;
				}
				else if ((nb[0] == 1 && nb[1] == 3)
					|| (nb[0] == 3 && nb[1] == 1))
				{
					return 4;
				}
				else if ((nb[0] == 2 && nb[1] == 3)
					|| (nb[0] == 3 && nb[1] == 2))
				{
					return 5;
				}
			}
		//}
	}

	Summit* GetBarycentreFromSummits(std::map<Summit*, float> m)
	{
		float x = 0, y = 0, z = 0;
		std::map<Summit*, float>::iterator it;
		for (it = m.begin(); it != m.end(); ++it)
		{
			x += it->first->getPoint().x * it->second;
			y += it->first->getPoint().y * it->second;
			z += it->first->getPoint().z * it->second;
		}
		Point p = Point(x, y, z);
		Summit *s = new Summit(p);
		return s;
	}

	Summit* GetBoundaryBarycentre(Edge* e)
	{
		std::map<Summit*, float> bary_summits;
		bary_summits[e->getSummitsConnected()->at(0)] = 9.0f / 16.0f;
		bary_summits[e->getSummitsConnected()->at(1)] = 9.0f / 16.0f;

		for (int j = 0; j < e->getSummitsConnected()->size(); ++j)		//for each summit of edge
		{
			Edge* edge = nullptr;
			for (int k = 0; k < e->getSummitsConnected()->at(j)->getEdgesConnected()->size(); ++k)		//for each adjacent face of summit
			{
				if (e->getSummitsConnected()->at(j)->getEdgesConnected()->at(k) != e
					&& e->getSummitsConnected()->at(j)->getEdgesConnected()->at(k)->getFacesConnected()->size() == 1)
				{
					edge = e->getSummitsConnected()->at(j)->getEdgesConnected()->at(k);
					break;
				}
			}
			bary_summits[GetOppositeSummit(e->getSummitsConnected()->at(j), edge)] = -1.0f / 16.0f;
		}

		return GetBarycentreFromSummits(bary_summits);
	}

	Summit* GetFullSchemeBarycentre(Edge* e)
	{
		std::map<Summit*, float> bary_summits;
		Summit* new_s;
		new_s = new Summit(e->getSummitsConnected()->at(0)->getPoint());
		bary_summits[new_s] = 1.0f / 2.0f;
		new_s = new Summit(e->getSummitsConnected()->at(1)->getPoint());
		bary_summits[new_s] = 1.0f / 2.0f;

		for (int i = 0; i < e->getFacesConnected()->size(); ++i)
		{
			for (int j = 0; j < e->getFacesConnected()->at(i)->getEdgesConnected()->size(); ++j)
			{
				if (e->getFacesConnected()->at(i)->getEdgesConnected()->at(j) == e)
				{
					continue;
				}

				for (int k = 0; k < e->getFacesConnected()->at(i)->getEdgesConnected()->at(j)->getFacesConnected()->size(); ++k)
				{
					if (e->getFacesConnected()->at(i)->getEdgesConnected()->at(j)->getFacesConnected()->at(k) == e->getFacesConnected()->at(i))
					{
						continue;
					}
					new_s = new Summit(GetOppositeSummit(e->getFacesConnected()->at(i)->getEdgesConnected()->at(j), e->getFacesConnected()->at(i)->getEdgesConnected()->at(j)->getFacesConnected()->at(k))->getPoint());
					bary_summits[new_s] = -1.0f / 16.0f;
				}
			}

			new_s = new Summit(GetOppositeSummit(e, e->getFacesConnected()->at(i))->getPoint());
			bary_summits[new_s] = 1.0f / 8.0f;
		}

		return GetBarycentreFromSummits(bary_summits);
	}

	Summit* GetSemiBoundaryBarycentre(Edge* e)
	{
		std::map<Summit*, float> bary_summits;
		bary_summits[e->getSummitsConnected()->at(0)] = 1.0f / 2.0f;
		bary_summits[e->getSummitsConnected()->at(1)] = 1.0f / 2.0f;

		return GetBarycentreFromSummits(bary_summits);
	}

	Summit* GetPyramidalBarycentre(Edge* e)
	{
		std::map<Summit*, float> bary_summits;
		bary_summits[e->getSummitsConnected()->at(0)] = 1.0f / 2.0f;
		bary_summits[e->getSummitsConnected()->at(1)] = 1.0f / 2.0f;

		for (int i = 0; i < e->getFacesConnected()->size(); ++i)
		{
			bool empty = false;
			for (int j = 0; j < e->getFacesConnected()->at(i)->getEdgesConnected()->size(); ++j)
			{
				if (e->getFacesConnected()->at(i)->getEdgesConnected()->at(j) != e)
				{
					continue;
				}
				else if (e->getFacesConnected()->at(i)->getEdgesConnected()->at(j)->getFacesConnected()->size() == 1)
				{
					empty = true;
					break;
				}

				for (int k = 0; k < e->getFacesConnected()->at(i)->getEdgesConnected()->at(j)->getFacesConnected()->size(); ++k)
				{
					if (e->getFacesConnected()->at(i)->getEdgesConnected()->at(j)->getFacesConnected()->at(k) == e->getFacesConnected()->at(i))
					{
						continue;
					}

					bary_summits[GetOppositeSummit(e->getFacesConnected()->at(i)->getEdgesConnected()->at(j), e->getFacesConnected()->at(i)->getEdgesConnected()->at(j)->getFacesConnected()->at(k))] = -1.0f / 8.0f;
				}
			}

			if (!empty)
			{
				bary_summits[GetOppositeSummit(e, e->getFacesConnected()->at(i))] = 1.0f / 4.0f;
			}
		}

		return GetBarycentreFromSummits(bary_summits);
	}

	Summit* GetUnbalancedBarycentre(Edge* e)
	{
		std::map<Summit*, float> bary_summits;

		Summit* common_summit = nullptr;
		int full_face = 0;
		for (int i = 0; i < e->getFacesConnected()->size(); ++i)
		{
			bool full = true;
			for (int j = 0; j < e->getFacesConnected()->at(i)->getEdgesConnected()->size(); ++j)
			{
				if (e->getFacesConnected()->at(i)->getEdgesConnected()->at(j)->getFacesConnected()->size() == 1)
				{
					full = false;
					break;
				}
			}

			if (!full)
			{
				for (int j = 0; j < e->getFacesConnected()->at(i)->getEdgesConnected()->size(); ++j)
				{
					if (e->getFacesConnected()->at(i)->getEdgesConnected()->at(j) == e
						|| e->getFacesConnected()->at(i)->getEdgesConnected()->at(j)->getFacesConnected()->size() == 1)
					{
						continue;
					}

					for (int k = 0; k < e->getFacesConnected()->at(i)->getEdgesConnected()->at(j)->getFacesConnected()->size(); ++k)
					{
						if (e->getFacesConnected()->at(i)->getEdgesConnected()->at(j)->getFacesConnected()->at(k) == e->getFacesConnected()->at(i))
						{
							continue;
						}

						bary_summits[GetOppositeSummit(e->getFacesConnected()->at(i)->getEdgesConnected()->at(j),
							e->getFacesConnected()->at(i)->getEdgesConnected()->at(j)->getFacesConnected()->at(k))] = -1 / 16;

						if (IsSummitInFace(e->getSummitsConnected()->at(0), e->getFacesConnected()->at(i)->getEdgesConnected()->at(j)->getFacesConnected()->at(k)))
						{
							common_summit = e->getSummitsConnected()->at(0);
							bary_summits[e->getSummitsConnected()->at(0)] = 5.0f / 8.0f;
							bary_summits[e->getSummitsConnected()->at(1)] = 3.0f / 8.0f;
						}
						else
						{
							common_summit = e->getSummitsConnected()->at(1);
							bary_summits[e->getSummitsConnected()->at(0)] = 3.0f / 8.0f;
							bary_summits[e->getSummitsConnected()->at(1)] = 5.0f / 8.0f;
						}
					}
				}

				bary_summits[GetOppositeSummit(e, e->getFacesConnected()->at(i))] = 1.0f / 16.0f;
			}
			else
			{
				full_face = i;
				continue;
			}
		}

		for (int j = 0; j < e->getFacesConnected()->at(full_face)->getEdgesConnected()->size(); ++j)
		{
			if (e->getFacesConnected()->at(full_face)->getEdgesConnected()->at(j) == e)
			{
				continue;
			}

			for (int k = 0; k < e->getFacesConnected()->at(full_face)->getEdgesConnected()->at(j)->getFacesConnected()->size(); ++k)
			{
				if (e->getFacesConnected()->at(full_face)->getEdgesConnected()->at(j)->getFacesConnected()->at(k) == e->getFacesConnected()->at(full_face))
				{
					continue;
				}

				if (IsSummitInFace(common_summit, e->getFacesConnected()->at(full_face)->getEdgesConnected()->at(j)->getFacesConnected()->at(k)))
				{
					bary_summits[GetOppositeSummit(e->getFacesConnected()->at(full_face)->getEdgesConnected()->at(j),
						e->getFacesConnected()->at(full_face)->getEdgesConnected()->at(j)->getFacesConnected()->at(k))] = -1.0f / 8.0f;
				}
				else
				{
					bary_summits[GetOppositeSummit(e->getFacesConnected()->at(full_face)->getEdgesConnected()->at(j),
						e->getFacesConnected()->at(full_face)->getEdgesConnected()->at(j)->getFacesConnected()->at(k))] = -1.0f / 16.0f;
				}
			}

			bary_summits[GetOppositeSummit(e, e->getFacesConnected()->at(full_face))] = 3.0f / 16.0f;
		}

		return GetBarycentreFromSummits(bary_summits);
	}

	bool IsSummitInFace(Summit* s, Face* f)
	{
		for (int i = 0; i < f->getSummitsConnected()->size(); ++i)
		{
			if (f->getSummitsConnected()->at(i) == s)
			{
				return true;
			}
		}

		return false;
	}

	Edge* EdgeExists(Edge* e, std::vector<Edge*>* edges)
	{
		for (int i = 0; i < edges->size(); ++i)
		{
			if ((e->getSummitsConnected()->at(0) == edges->at(i)->getSummitsConnected()->at(0) && e->getSummitsConnected()->at(1) == edges->at(i)->getSummitsConnected()->at(1))
				|| (e->getSummitsConnected()->at(1) == edges->at(i)->getSummitsConnected()->at(0) && e->getSummitsConnected()->at(0) == edges->at(i)->getSummitsConnected()->at(1)))
			{
				return edges->at(i);
			}
		}

		return nullptr;
	}

	Edge* FindEdgeFromSummits(Summit* s1, Summit* s2, std::vector<Edge*>* edges)
	{
		for (int i = 0; i < edges->size(); ++i)
		{
			if ((edges->at(i)->getSummitsConnected()->at(0) == s1 && edges->at(i)->getSummitsConnected()->at(1) == s2)
				|| (edges->at(i)->getSummitsConnected()->at(1) == s1 && edges->at(i)->getSummitsConnected()->at(0) == s2))
			{
				return edges->at(i);
			}
		}

		return nullptr;
	}

	Face* FindFaceFromSummits(Summit* s1, Summit* s2, Summit* s3, std::vector<Face*>* faces)
	{
		for (int i = 0; i < faces->size(); ++i)
		{
			if (faces->at(i)->getSummitsConnected()->at(0) == s1 || faces->at(i)->getSummitsConnected()->at(0) == s2 || faces->at(i)->getSummitsConnected()->at(0) == s3)
			{
				if (faces->at(i)->getSummitsConnected()->at(1) == s1 || faces->at(i)->getSummitsConnected()->at(1) == s2 || faces->at(i)->getSummitsConnected()->at(1) == s3)
				{
					if (faces->at(i)->getSummitsConnected()->at(2) == s1 || faces->at(i)->getSummitsConnected()->at(2) == s2 || faces->at(i)->getSummitsConnected()->at(2) == s3)
					{
						return faces->at(i);
					}
				}
			}
		}

		return nullptr;
	}
};
