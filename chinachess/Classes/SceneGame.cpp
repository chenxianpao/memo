#include "SceneGame.h"

SceneGame::SceneGame()
{
}


CCScene* SceneGame::scene(bool red)
{
    CCScene* scene = CCScene::create();

    SceneGame* layer = SceneGame::create(red);

    scene->addChild(layer);

    return scene;
}

//�Զ���create����
SceneGame* SceneGame::create(bool red)
{
    SceneGame* pRet = new SceneGame();
    
    if(pRet)
    {
        if(pRet->init(red))
        {
            pRet->autorelease();
        }
        else
        {
            pRet->release();
            return NULL;
        }
    }
    else
    {
        return NULL;
    }

    return pRet;
}

bool SceneGame::init(bool red)
{
    //���ø���CCLayer
    if(!CCLayer::init())
    {
        return false;
    }

     CCSize winSize = CCDirector::sharedDirector()->getWinSize();

     //�������̵�ƫ��ֵ
     this->_plateOffset = ccp(20,10);

     //�������ӵ�ƫ��ֵ
      this->_stoneOffset = ccp(60, 33);

      //�������ӵ�ֱ��Ϊ46
	  this->_d = 46;

      //��ʼ��ʱ��û��ѡ������
      _selectid = -1;

      //����ʱ���л����ӵ���ɫ
      _redTrun = true;

      //��������Ұں���
      _redSide = red;

     //��������
     CCSprite* desk = CCSprite::create("floor.jpg");
     this->addChild(desk);
     
     //�������ӵ�λ��
     desk->setPosition(ccp(winSize.width / 2, winSize.height / 2));
     
     //ѹ������
     desk->setScaleX(winSize.width / desk->getContentSize().width);
     desk->setScaleY(winSize.height / desk->getContentSize().height);
    
    
    //��������
    CCSprite* plate = CCSprite::create("background.png");
    this->addChild(plate);

    //�������Ϊ(0,0)
    plate->setAnchorPoint(CCPointZero);

    //�������̵�λ��
    plate->setPosition(_plateOffset);

    //ѹ�����̣�(���ڵĸ߶� - ƫ�Ƶ�y���� * 2) / ͼƬ�ĸ߶�
    plate->setScale((winSize.height -_plateOffset.y *2)/ plate->getContentSize().height);


    //������
    for(int i=0; i<32; i++)
    {
        //��������
        _s[i] = Stone::create(i, red);
        addChild(_s[i]);

        //�������ӵĳ�ʼλ��Ϊ���λ��
      _s[i]->setPosition(ccp(CCRANDOM_0_1() * winSize.width,
                             CCRANDOM_0_1() * winSize.height));

      //��������
      _s[i]->setVisible(false);
    }

    //����һ��ѡ���
    //��ѡ��ĳ�����ӵ�ʱ��,ѡ��������ѡ�õ�������
    _selectSprite = CCSprite::create("selected.png");
    addChild(_selectSprite);

    //����ѡ���
    _selectSprite->setVisible(false);

    //����ѡ�������ȼ�
    _selectSprite->setZOrder(1000);

    //ѹ��ѡ���
    _selectSprite->setScale(.8f);


    //����Menu
   CCMenu* menu = CCMenu::create();
   this->addChild(menu);

    //������ʼ��ť
    CCMenuItem* itemStart = CCMenuItemImage::create("start.jpg", "start.jpg",
                                           this, menu_selector(SceneGame::Start));
    menu->addChild(itemStart);
    itemStart->setPositionX(190);
    itemStart->setPositionY(120);

    //�����¾ְ�ť
    CCMenuItem* itemNew = CCMenuItemImage::create("new.jpg", "new.jpg",
                                           this, menu_selector(SceneGame::New));
    menu->addChild(itemNew);
    itemNew->setPositionX(itemStart->getPositionX());
    itemNew->setPositionY(itemStart->getPositionY() + 60);


   //�������尴ť
    CCMenuItem* item = CCMenuItemImage::create("regret.jpg", "regret.jpg",
                                               this, menu_selector(SceneGame::Back));
    menu->addChild(item);
    item->setPositionX(itemStart->getPositionX());
    item->setPositionY(itemStart->getPositionY() - 60);

    //������ͣ��ť
    CCMenuItem* itemPause = CCMenuItemImage::create("pause.jpg", "pause.jpg",
                                           this, menu_selector(SceneGame::Pause));
    menu->addChild(itemPause);
    itemPause->setPositionX(itemStart->getPositionX());
    itemPause->setPositionY(itemStart->getPositionY() - 60 - 60);

     //�����ѶȰ�ť
    CCMenuItem* itemDifficulty = CCMenuItemImage::create("difficulty.jpg", "difficulty.jpg",
                                           this, menu_selector(SceneGame::Difficulty));
    menu->addChild(itemDifficulty);
    itemDifficulty->setPositionX(itemStart->getPositionX());
    itemDifficulty->setPositionY(itemStart->getPositionY() - 60 - 60 - 60);

     //�������ű������ְ�ť
    CCMenuItem* itemVoice = CCMenuItemImage::create("openVolice.png", "closeVolice.png",
                                           this, menu_selector(SceneGame::Voice));
    menu->addChild(itemVoice);
    itemVoice->setPositionX(itemStart->getPositionX());
    itemVoice->setPositionY(itemStart->getPositionY() - 60 - 60 - 60 - 60);

    CCLog("x=%lf", itemStart->getPositionX());
    CCLog("y=%lf", itemStart->getPositionY() - 240);

    //����һ������,�����Դ�Ϊ���������  
    CCLabelTTF* label = CCLabelTTF::create("Voice", "Arial", 25);  
    addChild(label);  
 
   //�������ֵ�λ��  
   label->setPosition(ccp(winSize.width/2 + 120, winSize.height/2 - 120));  

   //�������ֵ���ɫ
   label->setColor(ccc3(0, 0, 0));


    //��������
    _steps = CCArray::create();
    _steps->retain();

    //����������ʾ��Ϸ���
    sprite = CCSprite::create("yingjiemian.png");
    addChild(sprite);
    sprite->setPosition(ccp(winSize.width / 2, winSize.height / 2));

    //���ؽ��
    sprite->setVisible(false);
    visible = false;

    //����
    setTouchEnabled(true);
    setTouchMode(kCCTouchesOneByOne);

    return true;
}

//ͨ�����ѡ�����ӣ�������
bool SceneGame::ccTouchBegan(CCTouch* pTouch, CCEvent* pEvent)
{
    CCObject* obj = (CCObject*)pTouch;

    //��ȡ������Ĵ�������
    CCPoint ptInWin = pTouch->getLocation();

  
     if(sprite->boundingBox().containsPoint(ptInWin) && visible == true)  
     {  
         //������Ϸ���
         sprite->setVisible(false);

         visible = false;

         New(obj);
     }

    int x, y;//���津�������������

    //ͨ��������Ĵ��������ȡ���̵�x�����y����
    if(!getClickPos(ptInWin, x, y))
    {
        return false;
    }

    //ͨ���������������е������ȡѡ�е����ӵ�id
    int clickid = getStone(x, y);
    //���������λ���������ӵ�ʱ��,clickidΪѡ�е����ӵ�id,��ʾ�����ѡ��
    //���������λ����û�����ӵ�ʱ��,clickidΪ-1,��ʾ���������

    //-1 == _selectid��ʾû��ѡ������
    if(-1 == _selectid)
    {
        setSelectId(clickid);
    }
    else
    {
        //�ƶ�����
        //��һ���������ƶ������ӵ�id
        //�ڶ���������ͨ���������λ���жϴ��������Ƿ�������
        //�������������������x����
        //���ĸ��������������y����
        //moveStoneִ������������ѡ�������
        //ѡ���ӣ���_selectid == clickidʱ����ʾѡ����idΪ_selectid������
        //�����ӣ���selectid != clickidʱ�� ��ʾ��idΪ_selectid�������ƶ���(x,y)���ڵ�λ����
        moveStone(_selectid, clickid, x, y);
    }

    // CCLog("_selectid=%d, clickid=%d", _selectid, clickid);
     //CCLog("x=%d, y=%d", x, y);

    return true;
}

//�õ�������������ϵ������
//���������λ���������ⷵ��false
//ͨ��������������������
bool SceneGame::getClickPos(CCPoint ptInWin, int &x, int &y)
{
    for(x=0; x<=8; x++)
    {
        for(y=0; y<=9; y++)
        {
            //���������ϵĸ����ڴ����ϵ�λ��
            CCPoint ptInPlate = getStonePos(x, y);

           // CCLog("ptInPlate.x=%lf   ptInPlate.y=%lf", ptInPlate.x,  ptInPlate.y);

            //Ѱ�����������λ�þ���С�����ӵİ뾶�ĸ���
            //����ҵ���,return true,���򷵻� return false
            if(ptInWin.getDistance(ptInPlate) < _d / 2)
            {
                return true;
            }
        }
    }

    return false;
}

//ͨ��������±��ȡ���ӵ�ID
//���������û������,����-1
int SceneGame::getStone(int x, int y)
{
    Stone* s;

    //����32������
    for(int i=0; i<32; i++)
    {
        s = _s[i];

        if(s->getX() == x && s->getY() == y && !s->getDead())
        {
            //�õ����ӵ�ID
            return s->getID();
        }
    }

	return -1;
}

//����ѡ���
void SceneGame::setSelectId(int id)
{
    if(-1 == id)
    {
        return;
    }

    //���û��ѡ�к���
    if(_s[id]->getRed() != _redTrun)
    {
        return;
    }

    //_selectidΪѡ�е����ӵ�id
    _selectid = id;

    //��ʾѡ���
    _selectSprite->setVisible(true);

    //ѡ����ڱ�ѡ�е���������ʾ
    _selectSprite->setPosition(_s[_selectid]->getPosition());
}


//�ƶ�����
//��һ���������ƶ������ӵ�id
//�ڶ���������ͨ���������λ���жϴ��������Ƿ�������
//�������������������x����
//���ĸ��������������y����
void SceneGame::moveStone(int moveId, int killId, int x, int y)
{
    //killId != -1��ʾ�������λ������һ������
    //_s[moveId]->getRed() == _s[killId]->getRed()��ʾ��������
    //�����Ӻ���������ӵ���ɫ��ͬ
    if(killId != -1 && _s[moveId]->getRed() == _s[killId]->getRed())
    {
        //����ѡ���
        setSelectId(killId);

        return;
    }

    //CCLog("killId=%d, moveId=%d", killId, moveId);
    //CCLog("_s[moveId]->getRed()=%d", _s[moveId]->getRed());
    

    //�������
   bool bCanMove =  canMove(moveId, killId, x, y);

   //���bCanMoveΪfalse
   //��������
   if(false == bCanMove)
   {
       return;
   }

    //����֮ǰ��¼���ӵ���Ϣ
   //��һ����������Ҫ�ƶ������ӵ�id
   //�ڶ���������ͨ���������λ���жϴ��������Ƿ�������
   //���������������ӵ�ǰ��λ�õ�x����
   //���ĸ����������ӵ�ǰ��λ�õ�y����
   //����������������ƶ����λ�õ�x����
   //�����������������ƶ����λ�õ�y����
    Step* step = Step::create(moveId, killId, _s[moveId]->getX(), _s[moveId]->getY(), x, y);
    
    //�����ӵ���Ϣ��ӵ�������
    _steps->addObject(step);

    //�������ӵ�����(�ƶ�����)
    _s[moveId]->setX(x);
    _s[moveId]->setY(y);

    //_s[moveId]->setPosition(getStonePos(x,y));
    //SetRealPos(_s[moveId]);

    //�����ƶ�����ʱ�Ķ���
    CCMoveTo* move = CCMoveTo::create(.5f, getStonePos(x, y));
    
    //�����ص�
    CCCallFuncND* call = CCCallFuncND::create(this, 
                            callfuncND_selector(SceneGame::moveComplete), 
                         (void*)(intptr_t)killId);

    //���ö�����ִ��˳��(���ƶ�����,����ûص�����)
    CCSequence* seq = CCSequence::create(move, call, NULL);
    
    //�����ƶ������ӵ����ȼ�
    _s[moveId]->setZOrder(_s[moveId]->getZOrder() + 1);

    //ִ�������ƶ��Ķ���
    _s[moveId]->runAction(seq);
}

//����������ת���ɴ�������
CCPoint SceneGame::getStonePos(int x, int y)
{
    int xx = x * _d;
    int yy = y * _d;

   return ccp(xx, yy) + _stoneOffset;
}


//ʵ�ֻ���
void SceneGame::Back(CCObject*)
{
    //�������е�Ԫ�ظ���Ϊ0ʱ
    //û����
    if(0 == _steps->count())
    {
        return;
    }

    //��ȡ�����е����һ��Ԫ��
    //��ȡ����ʱ�����һ�����ӵ���Ϣ
    Step* step = (Step*)_steps->lastObject();

   // �ָ����ӵ���Ϣ
    //������������ǰ��λ��x����
    _s[step->_moveid]->setX(step->_xFrom);

    //������������ǰ��λ��y����
    _s[step->_moveid]->setY(step->_yFrom);
    _s[step->_moveid]->setPosition(getStonePos(step->_xFrom, step->_yFrom));

    //�ָ��Ե�������
    if(step->_killid != -1)
    {
        //��ʾ�Ե�������
        _s[step->_killid]->setVisible(true);

        //����Ե�������
         _s[step->_killid]->setDead(false);
    }

    //�ƶ���һ�����
    //�л��ƶ������ӵ���ɫ
    _redTrun = ! _redTrun;

    //ɾ�������е����һ��Ԫ��
    //ɾ������ʱ���һ�����ӵ���Ϣ
    _steps->removeLastObject();
}


 //ʵ�ֿ�ʼ
void SceneGame::Start(CCObject*)
{
   

    //������
    for(int i=0; i<32; i++)
    {
        //��ʾ����
      _s[i]->setVisible(true);

      //�������ƶ���������ָ����λ��
      CCMoveTo* move = CCMoveTo::create(1, this->getStonePos(_s[i]->getX(), _s[i]->getY()));
      _s[i]->runAction(move);
    }
}


//ʵ���¾�
void SceneGame::New(CCObject* obj)
{
  //�õ����ڵĴ�С
    CCSize winSize = CCDirector::sharedDirector()->getWinSize();

     //���˶��ٲ���ͻڶ��ٲ�����
    for(int i = _steps->count(); i>0; i--)
    {
        Back(obj);
    }
    
     for(int i=0; i<32; i++)
    {
         //�������ӵĳ�ʼλ��Ϊ���λ��
        _s[i]->setPosition(ccp(CCRANDOM_0_1() * winSize.width,
                             CCRANDOM_0_1() * winSize.height));
      //��������
      _s[i]->setVisible(false); 
    }
}

//ʵ����ͣ
void SceneGame::Pause(CCObject*)
{
}


//ʵ���Ѷ�
void SceneGame::Difficulty(CCObject*)
{
}


//���ű�������
void SceneGame::Voice(CCObject*)
{
    static int i = 0;

    if(0 == i)
    {
        //���ű�������
        CocosDenshion::SimpleAudioEngine::sharedEngine()->playBackgroundMusic("floor.wav",true);
        
        i++;
    }
    else
    {
        //ֹͣ���ű�������
        CocosDenshion::SimpleAudioEngine::sharedEngine()->stopBackgroundMusic();

        i--;
    }
}



void SceneGame::moveComplete(CCNode* movestone, void* _killid)
{
    //�õ����ڵĴ�С
    CCSize winSize = CCDirector::sharedDirector()->getWinSize();

    //�������ȼ�
    movestone->setZOrder(movestone->getZOrder() - 1);

    int killid =  (int)(intptr_t)_killid;

    //�����������������
      if(-1 != killid)
    {
        //ɱ���������ϵ�����
        _s[killid]->setDead(true);

        //����ɱ��������
        _s[killid]->setVisible(false);

        //��ɱ������ʱ��,���¿�ʼ
        if(Stone::JIANG  == _s[killid]->getType())
        {
            //��ʾ��Ϸ���
            sprite->setVisible(true);

            //�������ȼ�
            sprite->setZOrder(1000);

            visible = true;
        }
    }
      
      //û��ѡ������
      _selectid = -1;

    //����ѡ���
    _selectSprite->setVisible(false);

    //�ƶ���һ�����
    //�л��ƶ������ӵ���ɫ
    _redTrun = ! _redTrun;
}


//�������
bool SceneGame::canMove(int moveid, int killid, int x, int y)
{
    //���ѡ�е�����
    Stone* s = _s[moveid];

    //���ӵ�����
    switch(s->getType())
    {
        //�����������
        case Stone::JIANG:
        {
            return canMoveJiang(moveid, killid, x, y);
        }
        break;

        //ʿ���������
        case Stone::SHI:
        {
            return canMoveShi(moveid, x, y);
        }
        break;

        //����������
        case Stone::XIANG:
        {
            return canMoveXiang(moveid, x, y);
        }
        break;
       
        //�����������
        case Stone::CHE:
        {
            return canMoveChe(moveid, x, y);
        }
        break;
       
        //����������
        case Stone::MA:
        {
            return canMoveMa(moveid, x, y);
        }
        break;
    
        //�ڵ��������
        case Stone::PAO:
        {
            return canMovePao(moveid, killid, x, y);
        }
        break;
     
        //�����������
        case Stone::BING:
        {
            return canMoveBing(moveid, x, y);
        }
        break;

        default:
        {
            break;
        }
    }

    return false;
}


//�����������
bool SceneGame::canMoveJiang(int moveid, int killid, int x, int y)
{
   Stone* skill = _s[killid];

      //�����������
    //1��һ����һ��
    //2�����ܳ��Ź���


      CCLog("x=%d, y=%d", x, y);
     CCLog("moveid=%d, killid=%d", moveid, killid);

    //���Ķ�ɱ
    if(skill->getType() == Stone::JIANG)
    {
        return canMoveChe(moveid, x, y);
    }

    //ͨ�����ӵ�ID�õ�����
    Stone* s = _s[moveid];

    //��ý���ǰ��λ��
    int xo = s->getX();
    int yo = s->getY();

    //��ý��ߵĸ���
    //(x,y)��ʾ���ߵ���λ��
    int xoff = abs(xo - x);
    int yoff = abs(yo - y);
    
    int d = xoff*10 + yoff;

    //�߽���ʱ�����������
    //xoff=1, yoff=0�������������
    //xoff=0, yoff=1������ǰ�����
    if(d != 1 && d != 10)
    {
        return false;
    }

    //�жϽ��Ƿ���˾Ź�
    //��ɫ�Ľ��ͺ�ɫ�Ľ���x����ķ�Χ����3<=x<=5
    if(x<3 || x>5)
    {
        return false;
    }

    //�����ҵ������Ǻ���
    if(_redSide == s->getRed())
    {
        //�жϽ��Ƿ���˾Ź�
        if(y<0 || y>2)
        {
            return false;
        }
    }
    else//�жϺ�ɫ�Ľ��ķ�Χ
    {
        //�жϽ��Ƿ���˾Ź�
        if(y>9 || y<7)
        {
            return false;
        }
    }

    return true;
}


//ʿ���������
bool SceneGame::canMoveShi(int moveid, int x, int y)
{
    //ʿ���������
    //1��һ����һ��
    //2�����ܳ��Ź���
   //3��б����

     //ͨ�����ӵ�ID�õ�����
    Stone* s = _s[moveid];

    //���������ǰ��λ��
    int xo = s->getX();
    int yo = s->getY();

    //������ߵĸ���
    //(x,y)��ʾ���ߵ���λ��
    int xoff = abs(xo - x);
    int yoff = abs(yo - y);

    int d = xoff*10 + yoff;

    //ʿÿ��һ��x������1��,y������1��
    //���ߵĸ�������1��ʱ
    //����false
    if(d != 11)
    {
        return false;
    }

     //�ж�ʿ�Ƿ���˾Ź�
    //��ɫ��ʿ�ͺ�ɫ��ʿ��x����ķ�Χ����3<=x<=5
    if(x<3 || x>5)
    {
        return false;
    }

    //�����ҵ������Ǻ���
    if(_redSide == s->getRed())
    {
        //�ж�ʿ�Ƿ���˾Ź�
        if(y<0 || y>2)
        {
            return false;
        }
    }
    else//�жϺ�ɫ��ʿ�ķ�Χ
    {
        //�ж�ʿ�Ƿ���˾Ź�
        if(y>9 || y<7)
        {
            return false;
        }
    }

    return true;
}


//����������
bool SceneGame::canMoveXiang(int moveid, int x, int y)
{
     //����������
    //ÿ��һ��x�ƶ�2��,y�ƶ�2��
    //���ܹ���


    //ͨ�����ӵ�ID�õ�����
    Stone* s = _s[moveid];

    //���������ǰ��λ��
    int xo = s->getX();
    int yo = s->getY();

    //������ߵĸ���
    //(x,y)��ʾ���ߵ���λ��
    int xoff = abs(xo - x);
    int yoff = abs(yo - y);

    int d = xoff*10 + yoff;

    //��ÿһ��x������2����,y������2��
    //���ߵĸ�������2��ʱ
    //����false
    if(d != 22)
    {
        return false;
    }

    //��������������е�����
    int xm = (xo + x) / 2;
    int ym = (yo + y) / 2;

    //�õ�(xm,ym)�ϵ�����
    int id = getStone(xm, ym);

    //��(xm,ym)�������ӵ�ʱ��
    if(id != -1)
    {
        //��������
        return false;
    }

      //�����಻�ܹ���
     //�����ҵ������Ǻ���
    if(_redSide == s->getRed())
    {
        //�ж����Ƿ���˺�
        if(y > 4)
        {
            return false;
        }
    }
    else//�жϺ�ɫ����ķ�Χ
    {
         //�ж����Ƿ���˺�
        if(y < 5)
        {
            return false;
        }
    }

    return true;
}


//�����������
bool SceneGame::canMoveChe(int moveid, int x, int y)
{
    //ͨ�����ӵ�ID�õ�����
    Stone* s = _s[moveid];

    //��ó�����ǰ��λ��
    int xo = s->getX();
    int yo = s->getY();

    //������֮�������ӵ�ʱ�򳵲�����
    if(getStoneCount(xo,yo,x,y) != 0)
    {
        return false;
    }

    return true;
}


//����������
bool SceneGame::canMoveMa(int moveid, int x, int y)
{
    //ͨ�����ӵ�ID�õ�����
    Stone* s = _s[moveid];

     //���������ǰ��λ��
    int xo = s->getX();
    int yo = s->getY();

    //CCLog("xo=%d", xo);
    //CCLog("yo=%d", yo);
    
     //������ߵĸ���
    //(x,y)��ʾ���ߵ���λ��
    //�������������
    //��һ�������������ǰ�������1�����������������2��
    //�ڶ�����������������������1��������ǰ�������2��
    int xoff = abs(xo-x);
    int yoff = abs(yo-y);

    //CCLog("x=%d", x);
    //CCLog("y=%d", y);
    
    int d = xoff*10 + yoff;

    //CCLog("d=%d", d);
    
    if(d != 12 && d != 21)     
    {
        return false;
    }

    int xm, ym;//��¼��ŵ�����
   
    if(d == 12)//�����ߵ��ǵ�һ�����
    {
        xm = xo;//��ŵ��x����Ϊ����ǰ���x����
        ym = (yo + y) / 2;//��ŵ��y����Ϊ����ǰ���y�������������y������е�����
    }
    else//�����ߵ��ǵڶ������
    {
        xm = (xo + x) / 2;//��ŵ��x����Ϊ����ǰ���x�������������x������е�����
        ym = yo;;//��ŵ��y����Ϊ����ǰ���y����
    }

    //CCLog("xm=%d", xm);
    //CCLog("ym=%d", ym);
    
    //����ŵ�������ʱ,������
    if(getStone(xm, ym) != -1) 
    {
        return false;
    }

    return true;
}


//�ڵ��������
bool SceneGame::canMovePao(int moveid, int killid, int x, int y)
{
    //ͨ�����ӵ�ID�õ�����
    Stone* s = _s[moveid];

    //���������ǰ��λ��
    int xo = s->getX();
    int yo = s->getY();

    //������������һ������
    //��������֮��ֻ��һ�����ӵ�ʱ��
    //�ڳԵ��������ϵ�����
    if(killid != -1 && this->getStoneCount(xo,yo,x,y) == 1)
    {
        return true;
    }

    if(killid == -1 && this->getStoneCount(xo, yo, x, y) == 0) 
    {
        return true;
    }

    return false;
}


//�����������
bool SceneGame::canMoveBing(int moveid, int x, int y)
{
     //�����������
    //1��һ����һ��
    //2��ǰ��һ����ܺ���
    //3�����Ӻ�ſ��������ƶ�

    //ͨ�����ӵ�ID�õ�����
    Stone* s = _s[moveid];

    //��ý���ǰ��λ��
    int xo = s->getX();
    int yo = s->getY();

    //��ñ��ߵĸ���
    //(x,y)��ʾ���ߵ���λ��
    int xoff = abs(xo - x);
    int yoff = abs(yo - y);
    
    int d = xoff*10 + yoff;

    //�߽���ʱ�����������
    //xoff=1, yoff=0�������������
    //xoff=0, yoff=1������ǰ�����
    if(d != 1 && d != 10)
    {
        return false;
    }

     //�����ҵ������Ǻ���
    if(_redSide == s->getRed())
    {
        //���ƺ�ɫ�ı����ܺ���
        if(y < yo)
        {
            return false;
        }

        //��ɫ�ı�û�й��Ӳ��������ƶ�
        if(yo <= 4 && y == yo)
        {
            return false;
        }
    }
    else//�жϺ�ɫ�ı�
    {
       //���ƺ�ɫ�ı����ܺ���
        if(y > yo)
        {
            return false;
        }

         //��ɫ�ı�û�й��Ӳ��������ƶ�
        if(yo >= 5 && y == yo)
        {
            return false;
        }
    }

    return true;
}


///����(xo,yo)��(x,y)֮���������
//���������Ϊ-1,��ʾ(xo,yo)��(x,y)����һ��ֱ����
int SceneGame::getStoneCount(int xo, int yo, int x, int y)
{
    int ret = 0;//��¼����֮������ӵĸ���

    //(xo,yo)��(x,y)����ͬһ��ֱ����
    if(xo != x && yo != y)
    {
        return -1;
    }

    //(xo,yo)��(x,y)��ͬһ����
    if(xo == x && yo == y)
    {
        return -1;
    }

    //������ͬһ��������
    if(xo == x)
    {
        //minΪ��������y������С�ĵ��y����
        int min = yo < y ? yo : y;

        //maxΪ��������y�������ĵ��y����
        int max = yo > y ? yo : y;

        //����ͬһ������������֮���������
        for(int yy=min+1; yy<max; yy++)
        {
            //������֮�������ӵ�ʱ��
            if(getStone(x,yy) != -1)
            {
                ++ret;//��������1
            }
        }
    }
    else//������ͬһ��������yo == y
    {
         //minΪ��������x������С�ĵ��x����
        int min = xo < x ? xo : x;

        //maxΪ��������x�������ĵ��x����
        int max = xo > x ? xo : x;

        //����ͬһ������������֮���������
        for(int xx=min+1; xx<max; xx++)
        {
             //������֮�������ӵ�ʱ��
            if(getStone(xx,y) != -1)
            {
                ++ret;//��������1
            }
        }
    }

    //��������֮���������
    return ret;
}
