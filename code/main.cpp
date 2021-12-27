#include <iostream>
#include "HoleFiller.h"
#include <spdlog/sinks/stdout_sinks.h>
//#include "MeshSmooth.hpp"
//#include "MeshSmooth.h"
#include "MeshDoctor.h"
#include "RBF_interpolation.hpp"
#include <SplineInterpolation.hpp>
#include "FindInnerBoundary.h"
using namespace std;
using namespace MeshFix;

void TransferJsonToVertical(std::string fileName);
//TransferJsonToVertical("C:/Users/st/Desktop/ToothRelated/AliSegmentData/data/lower_label.txt");
//1.spdlog 日志开源库
//2.OpenMesh 非结构网格开源库
//3.Hermite 插值我集成得到HoleFiller里面了
int main(int argc, char*argv[])
{
	Triangle_mesh inputMesh;
	bool inputFlag = OpenMesh::IO::read_mesh(inputMesh, "../../DataDir/lower.stl");//往上跳两级目录找到输入数据文件夹
	if (!inputFlag)
		return -1;
	FindInnerBoundary findBoundary(inputMesh);
	findBoundary.SetFaceLabelMap("../../DataDir/lower_label.txt");
	findBoundary.TraverseFaces();
	return system("pause");
}

//把json文件转换为垂直的
void TransferJsonToVertical(std::string fileName)
{
	ifstream ifile(fileName);
	//将文件读入到ostringstream对象buf中
	ostringstream buf;
	char ch;
	while (buf && ifile.get(ch))
	{
		switch (ch)
		{
		case ',':
			continue;
		case ' ':
			ch = '\n';
		default:
			buf.put(ch);;
		}
		
	}
	//返回与流对象buf关联的字符串
	//std::string a = buf.str();
	std::ofstream m_OutFile("D:/ttttttt.txt");
	m_OutFile << buf.str();
	m_OutFile.close();
	
}