#include "Stone.h"
 
Stone::Stone()
{
}

//��������
//��һ�����������ӵ�����
//�ڶ������������ӵ���ɫ
Stone* Stone::create(int id, bool red)
{
    Stone* s = new Stone();
    
    s->init(id, red);
        
    s->autorelease();
      
    return s;
}

//�ڷ�����
//�����������ж���Ϸ���ѡ�������Ϸ����ʱ�Ƿ�
//���˺�ɫ������
void Stone::reset(bool red)
{
    //�հ����ӵ�ʱ������û���Ե�
   this->setDead(false);
   
   if(red)//����ҵ��ɫ�����ӽ�����Ϸ����ʱ
    {//��ҵ����ӵ���ɫΪ��ɫ

         if(_id < 16)//�ڷź�ɫ������
        {
            //�������ӵ�x����
            this->setX(_initPos[_id]._x);

            //�������ӵ�y����
            this->setY(_initPos[_id]._y);
         }
        else//�ڷź�ɫ������
        {
            this->setX(8 - _initPos[_id-16]._x);

            //����������ͬ��ɫ��ͬ�����ӵ�y�������9
            //�磺��ɫ�ĳ��ͺ�ɫ�ĳ���y�������9
            this->setY(9 - _initPos[_id-16]._y);
        }
    }
    else//����ҵ��ɫ�����ӽ�����Ϸ����ʱ
    {//��ҵ����ӵ���ɫ�Ǻ�ɫ

          if(_id < 16)//�ڷź�ɫ������
        {
            this->setX(8 - _initPos[_id]._x);
            this->setY(9 - _initPos[_id]._y);
        }
        else//�ڷź�ɫ������
        {
            //����������ͬ��ɫ��ͬ�����ӵ�id���16
            //�磺��ɫ�ĳ��ͺ�ɫ�ĳ���id���16
            this->setX(_initPos[_id-16]._x);
            this->setY(_initPos[_id-16]._y);
        }
    }
}

//��ʼ������
//�����һ���Գ�ʼ��
bool Stone::init(int id, bool red)
{
     _id = id;//������ӵ�id

    //�����ӵ�idС��16ʱ,�����Ǻ�ɫ��
    _red = _id < 16;

    //��ʼ����ɫ������
     if(_id < 16)
     {
         _type = _initPos[_id]._type;
     }
     else//��ʼ����ɫ������
     {
         //����(��������)��ͬ��ɫ��ͬ���������ӵ�id���16
         _type = _initPos[_id-16]._type;
     }

    const char* stonePic[14] = {
            "rshuai.png", //(��ɫ)˧
            "rshi.png",   //(��ɫ)ʿ
            "rxiang.png", //(��ɫ)��
            "rche.png",   //(��ɫ)��
            "rma.png",    //(��ɫ)��
            "rpao.png",   //(��ɫ)��
            "rbing.png",  //(��ɫ)��

            "bjiang.png", //(��ɫ)��
            "bshi.png",   //(��ɫ)ʿ
            "bxiang.png", //(��ɫ)��
            "bche.png",   //(��ɫ)��
            "bma.png",    //(��ɫ)��
            "bpao.png",   //(��ɫ)��
            "bzu.png"     //(��ɫ)��
        };

    //����ͼƬ���±�
    //��������Ǻ�ɫ�� idx = _type
    //��������Ǻ�ɫ�� idx = 7 + _type
    //�������ӵ���ɫ���±���� 7
   int idx = (_red ? 0 : 1) * 7 + _type;

    //��������(��������)
    CCSprite::initWithFile(stonePic[idx]);

    //ѹ������
    setScale(.8f);

    //������(�������ӵ�λ��)
    reset(red);

    return true;
}

//����һ�ű�
Stone::InitPos Stone::_initPos[16] =
{
    //����λ��(0,0)
    {0, 0, Stone::CHE},

    //���λ��(1,0)
    {1, 0, Stone::MA},

    //���λ��(2,0)
    {2, 0, Stone::XIANG},

    //ʿ��λ��(3,0)
    {3, 0, Stone::SHI},

    //����λ��(4,0)
    {4, 0, Stone::JIANG},

     //ʿ��λ��(5,0)
    {5, 0, Stone::SHI},

    //���λ��(6,0)
    {6, 0, Stone::XIANG},

     //���λ��(7,0)
    {7, 0, Stone::MA},

     //����λ��(8,0)
    {8, 0, Stone::CHE},

    //�ڵ�λ��(1,2)
    {1, 2, Stone::PAO},

     //�ڵ�λ��(7,2)
    {7, 2, Stone::PAO},

    //����λ��(0,3)
    {0, 3, Stone::BING},

     //����λ��(2,3)
    {2, 3, Stone::BING},

     //����λ��(4,3)
    {4, 3, Stone::BING},

     //����λ��(6,3)
    {6, 3, Stone::BING},

     //����λ��(8,3)
    {8, 3, Stone::BING},
};