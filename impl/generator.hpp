#pragma once
#include <memory>
#include <vector>
#include <set>
#include "common.hpp"
#include "objectlib.hpp"
#include "monsterlib.hpp"

namespace snakegen
{

const std::size_t XL = 216;

enum class Presence
{
    Free = 0,
    Occupied
};

template<class T>
class Matrix
{
public:
    Matrix(size_t rows, size_t cols)
    : mRows(rows),
    mCols(cols),
    mData(rows * cols)
    {
    }

    T& operator()(size_t i, size_t j)
    {
        return mData[i * mCols + j];
    }

    T operator()(size_t i, size_t j) const
    {
        return mData[i * mCols + j];
    }

    int width() const
    {
        return mRows;
    }

    int heigth() const
    {
        return mRows;
    }

    bool contains(int x, int y) const
    {
        if(x < 0 && y < 0)
        {
            return false;
        }
        else
        {
            return static_cast<std::size_t>(x) < mRows && static_cast<std::size_t>(y) < mCols;
        }
    }

    void setAll(const T t)
    {
        std::fill(mData.begin(), mData.end(), t);
    }

    void merge(const Matrix<T>& other)
    {
        for(size_t i = 0; i < mData.size(); ++i)
        {
            mData[i] = static_cast<T>(std::max(mData[i], other.mData[i]));
        }
    }

private:
    size_t mRows;
    size_t mCols;
    std::vector<T> mData;
};

struct Zone
{
    Polygon polygon;
    std::vector<Coord> exitPoints;
    bool isStart = false;
};

struct Connection
{
    Connection(const std::size_t zoneA_, const std::size_t zoneB_)
    {
        if(zoneA_ > zoneB_)
        {
            zoneA = zoneB_;
            zoneB = zoneA_;
        }
        else
        {
            zoneA = zoneA_;
            zoneB = zoneB_;
        }
    }

    bool operator== (const Connection& other) const
    {
        return (zoneA == other.zoneA && zoneB == other.zoneB);
    }

    bool operator< (const Connection& other) const
    {
        return std::make_pair(zoneA, zoneB) < std::make_pair(other.zoneA, other.zoneB);
    }

    std::size_t first() const
    {
        return zoneA;
    }

    std::size_t second() const
    {
        return zoneB;
    }

private:
    std::size_t zoneA;
    std::size_t zoneB;
};

struct Map
{
    Map(const std::size_t size_)
    : size(size_)
    , obstacleLayer(size, size)
    , availabilityLayer(size, size)
    , roadLayer(size, size)
    , distanceLayer(size, size)
    {

    }

    int size = 0;
    Matrix<Presence> obstacleLayer;
    Matrix<Presence> availabilityLayer;
    Matrix<Presence> roadLayer;
    Matrix<double> distanceLayer;
    std::vector<Zone> zones;
    std::set<Connection> connections;
    std::vector<Object> objects;
    std::vector<Monster> monsters;
};

std::unique_ptr<Map> genMap(const std::size_t size, const std::size_t seed);

}