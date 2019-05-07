#include "SceneStart.h"
#include "SceneGame.h"

CCScene* SceneStart::scene()
{
    CCScene* scene = CCScene::create();
	if (!scene)
	{
		return NULL;
	}

    SceneStart* layer = SceneStart::create();
	if (!layer)
	{
		return NULL;
	}
    scene->addChild(layer);

    return scene;
}

bool SceneStart::init()
{
    CCLayer::init();

    CCSize winSize = CCDirector::sharedDirector()->getWinSize();

    CCSprite* bkr = CCSprite::create("bkg2.png");
    addChild(bkr);

    CCSprite* bkb = CCSprite::create("bkg1.png");
    addChild(bkb);

    bkr->setPosition(ccp(winSize.width / 2 - 100, winSize.height / 2));
    bkb->setPosition(ccp(winSize.width / 2 + 100, winSize.height / 2));

    _red = bkr;
    _black = bkb;

    //ע�ᴥ���¼�
    setTouchEnabled(true);
    setTouchMode(kCCTouchesOneByOne);
   
    return true;
}

bool SceneStart::ccTouchBegan(CCTouch* pTouch, CCEvent* pEvent)
{
    return true;
}

void SceneStart::ccTouchEnded(CCTouch* pTouch, CCEvent* pEvent)
{
     CCSize winSize = CCDirector::sharedDirector()->getWinSize();

    //��ô������λ��(����)
    CCPoint ptClick = pTouch->getLocation();

    //�����ж��Ƿ����������
    bool bClickStone = false;

    //�����к�ɫ�����ӵ�ʱ��(�������λ���ں�ɫ���������ڵķ�Χ��)
    if(_red->boundingBox().containsPoint(ptClick))
    {
        //�����˺�ɫ������
        this->_selected = true;

        //����������
        bClickStone = true;
    }
    //�����к�ɫ���ӵ�ʱ��(�������λ���ں�ɫ�������ڵķ�Χ��)
    else if(_black->boundingBox().containsPoint(ptClick))
    {
        //û���к�ɫ����
        this->_selected = false;

        //����������
        bClickStone = true;
    }

    //�����������ӵ�ʱ��
    if(bClickStone)
    {
        //�ƶ�����
        CCMoveTo* moveTo1 = CCMoveTo::create(1, ccp(winSize.width / 2, winSize.height / 2));
        CCMoveTo* moveTo2 = CCMoveTo::create(1, ccp(winSize.width / 2, winSize.height / 2));
    
        //��ת����
        CCRotateBy* rotate1 =  CCRotateBy::create(1, 360);
        CCRotateBy* rotate2 =  CCRotateBy::create(1, -360);

        //��ת���ƶ�ͬʱִ��
        CCSpawn* spawn1 = CCSpawn::create(moveTo1, rotate1, NULL);
        CCSpawn* spawn2 = CCSpawn::create(moveTo2, rotate2, NULL);

        //ִ���ж���
        _red->runAction(spawn1);
        _black->runAction(spawn2);

        //������ʱ��
         scheduleUpdate();
    }
}

void SceneStart::update(float)
{
    //��ȡ�������ӵ�x����
    float x1 = _red->getPositionX();
    float x2 = _black->getPositionX();

    //����ɫ�����Ӻͺ�ɫ��������ײ��
    //�������ӵľ���С�����ӵ�ֱ��
    //getContentSize().width������ӵĿ��(���ӵ�ֱ��)
    if(abs(x1 - x2) <= _red->getContentSize().width)
    {
        //������Ϸ
        CCDirector::sharedDirector()->replaceScene(SceneGame::scene(this->_selected));
    }
}