#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#include "Vector.h"
#include "ObjParser.h"
#include "VolumeGeometry.h"

namespace lux
{

const void ObjParser::split(const std::string& s, char delim, std::vector<std::string>& tokens)
{
	std::stringstream ss(s);
	std::string token;
	//std::cout << "Split " << s << " into ";
	while(std::getline(ss, token, delim))
	{
		if(!token.empty())
		{
			//std::cout<<token << " ";
			tokens.push_back(token);
		}
	}

	//std::cout <<std::endl;
}

void ObjParser::load(const std::string& filename, Mesh m)
{
	std::ifstream infile(filename);
	if(!infile.is_open())
	{
		std::cout << "File open failed" << std::endl;
		return;
	}

	std::string line;
	while(std::getline(infile, line))
	{
		if( IsVertex(line) )
		{
			Vector P;
			std::stringstream ss(line.substr(2));
			ss >> P[0] >> P[1] >> P[2];
			m->addVertex(P);
		}
		/*
		if( IsTexture(line) )
		{
			Vector P;
			std::stringstream ss(line.substr(2));
			ss >> P[0] >> P[1] >> P[2];
			m->addTextureCoordinate(P);
		}
		if( IsNormal(line) )
		{
			Vector P;
			std::stringstream ss(line.substr(2));
			ss >> P[0] >> P[1] >> P[2];
		}
		*/
		else if( IsFace(line) )
		{
			//std::cout << "line = " << line << std::endl;
			//	     {P0,t0,n0  P1,t1,n1  P2,t2,n2}
			int arr[9] = {-1,-1,-1, -1,-1,-1, -1,-1,-1};
			stringstream ss(line.substr(2));
			std::string vert;
			int index=0;
			while(ss >> vert)
			{
				std::vector<std::string> tokens;
				split(vert, '/', tokens);
				//FIX this mess later. Want to be able to keep texture/normal indices. Probable need to simplify
				/*
				std::cout << "size" <<tokens.size() << std::endl;
				for(size_t i=0; i<tokens.size() && i<3; ++i)
				{
					arr[index++] = tokens[i].empty() ? -1 : (std::stoi(tokens[i]) - 1);
				}*/
				arr[index] = std::stoi(tokens[0] )-1;
				index++;
			}
			
			//std::cout << "arr: ";
			//for(int i=0; i<9;i++) std::cout << arr[i] << " ";
			//std::cout << std::endl;
			m->addFace(arr[0], arr[1], arr[2]);
		}
	}
}

}
