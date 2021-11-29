#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

from sqlalchemy import Column, String, create_engine, MetaData, ForeignKey, Integer
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

Base = declarative_base()
meta = MetaData()


# 定义User对象
class Student(Base):
    __tablename__ = 'student'
    id = Column(String(20), primary_key=True)
    name = Column(String(20))


class Family(Base):
    __tablename__ = 'family'
    id = Column(String(20), primary_key=True)
    member = Column(Integer)
    student_id = Column(String(20), ForeignKey('student.id'))


class House(Base):
    __tablename__ = 'house'
    id = Column(String(20), primary_key=True)
    location = Column(String(100))
    family_id = Column(String(20), ForeignKey('family.id'))


class Car(Base):
    __tablename__ = 'car'
    id = Column(String(20), primary_key=True)
    name = Column(String(100))
    family_id = Column(String(20))


def create_fun():




    # 初始数据库连接
    #engine = create_engine('mysql+pymsql://huanggy:654321@localhost:3306/taxdb?charset=utf8', echo=True)
    HOST = 'localhost'
    PORT = 3306
    USERNAME = 'huanggy'
    PASSWORD = '654321'
    DB = 'taxdb'

    # dialect + driver://username:passwor@host:port/database
    #DB_URI = f'mysql+pymysql://{USERNAME}:{PASSWORD}@{HOST}:{PORT}/{DB}'
    #engine = create_engine(DB_URI, encoding="utf-8", echo=True)

    engine = create_engine(f"mysql+pymysql://{USERNAME}:{PASSWORD}@{HOST}:{PORT}/{DB}", encoding='utf-8',echo=True)

    # 创建DBsession
    DBSession = sessionmaker(bind=engine)

    # 创建session会话，数据库操作的基石。
    session = DBSession()

    # 在数据库中创建表user
    Student.metadata.create_all(bind=engine)
    Family.metadata.create_all(bind=engine)
    House.metadata.create_all(bind=engine)

    # 插入数据
    stu_one = Student(id='1', name='wukong')
    stu_two = Student(id='2', name='beijita')
    stu_three = Student(id='3', name='bike')
    #stu_four = Student(id='4', name='')

    # 提交数据到session

    session.add(stu_one)
    session.add(stu_two)
    session.add(stu_three)

    # session.add(stu_four)
    session.commit()

    family_one = Family(id='1', member=7, student_id='1')
    family_two = Family(id='2', member=5, student_id='2')
    family_three = Family(id='3', member=8, student_id='3')

    session.add(family_one)
    session.add(family_two)
    session.add(family_three)
    session.commit()

    house_one = House(id='1', location='earth', family_id='1')
    house_two = House(id='2', location='beijitaxing', family_id='2')
    house_three = House(id='3', location='meikexingren', family_id='3')
    house_four = House(id='4', location='earth', family_id='3')

    session.add(house_one)
    session.add(house_two)
    session.add(house_three)
    session.add(house_four)
    session.commit()

    car_one = Car(id='1', name='jindouyun', family_id='1')
    car_two = Car(id='2', name='benchi', family_id='2')
    car_three = Car(id='3', name='BWM', family_id='3')

    session.add(car_one)
    session.add(car_two)
    session.add(car_three)
    # 提交到数据库
    session.commit()

    session.close()


if __name__ == '__main__':
    #create_fun()
    #查询
    engine = create_engine(f"mysql+pymysql://huanggy:654321@localhost:3306/taxdb", encoding='utf-8',echo=True)

    # 创建DBsession
    DBSession = sessionmaker(bind=engine)

    # 创建session会话，数据库操作的基石。
    session = DBSession()
    result = session.query(Student).join(Family).filter(Family.member > 6)
    print('------->',result)

    result_two = session.query(Student).join(Family).join(House).filter(House.location == 'meikexingren')
    print('=======>',result_two)


