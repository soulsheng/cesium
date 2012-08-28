/*global defineSuite*/
defineSuite([
         'Core/BoundingRectangle',
         'Core/Cartesian2',
         'Core/Ellipsoid',
         'Core/EquidistantCylindricalProjection',
         'Core/Extent',
         'Core/Intersect'
     ], function(
         BoundingRectangle,
         Cartesian2,
         Ellipsoid,
         EquidistantCylindricalProjection,
         Extent,
         Intersect) {
    "use strict";
    /*global jasmine,describe,xdescribe,it,xit,expect,beforeEach,afterEach,beforeAll,afterAll,spyOn,runs,waits,waitsFor*/

    it('default constructor sets expected values', function() {
        var rectangle = new BoundingRectangle();
        expect(rectangle.x).toEqual(0.0);
        expect(rectangle.y).toEqual(0.0);
        expect(rectangle.width).toEqual(0.0);
        expect(rectangle.height).toEqual(0.0);
    });

    it('constructor sets expected parameters', function() {
        var rectangle = new BoundingRectangle(1.0, 2.0, 3.0, 4.0);
        expect(rectangle.x).toEqual(1.0);
        expect(rectangle.y).toEqual(2.0);
        expect(rectangle.width).toEqual(3.0);
        expect(rectangle.height).toEqual(4.0);
    });

    it('clone without a result parameter', function() {
        var rectangle = new BoundingRectangle(1.0, 2.0, 3.0, 4.0);
        var result = rectangle.clone();
        expect(rectangle).toNotBe(result);
        expect(rectangle).toEqual(result);
    });

    it('clone with a result parameter', function() {
        var rectangle = new BoundingRectangle(1.0, 2.0, 3.0, 4.0);
        var result = new BoundingRectangle(6.0, 7.0, 8.0, 9.0);
        var returnedResult = rectangle.clone(result);
        expect(result).toNotBe(rectangle);
        expect(result).toEqual(rectangle);
        expect(result).toBe(returnedResult);
    });

    it('clone works with "this" result parameter', function() {
        var rectangle = new BoundingRectangle(1.0, 2.0, 3.0, 4.0);
        var returnedResult = rectangle.clone(rectangle);
        expect(rectangle.x).toEqual(1.0);
        expect(rectangle.y).toEqual(2.0);
        expect(rectangle.width).toEqual(3.0);
        expect(rectangle.height).toEqual(4.0);
        expect(rectangle).toBe(returnedResult);
    });

    it('equals', function() {
        var rectangle = new BoundingRectangle(1.0, 2.0, 3.0, 4.0);
        expect(rectangle.equals(new BoundingRectangle(1.0, 2.0, 3.0, 4.0))).toEqual(true);
        expect(rectangle.equals(new BoundingRectangle(5.0, 2.0, 3.0, 4.0))).toEqual(false);
        expect(rectangle.equals(new BoundingRectangle(1.0, 6.0, 3.0, 4.0))).toEqual(false);
        expect(rectangle.equals(new BoundingRectangle(1.0, 2.0, 7.0, 4.0))).toEqual(false);
        expect(rectangle.equals(new BoundingRectangle(1.0, 2.0, 3.0, 8.0))).toEqual(false);
        expect(rectangle.equals(undefined)).toEqual(false);
    });

    var positions = [new Cartesian2(3, -1),
                     new Cartesian2(2, -2),
                     new Cartesian2(1, -3),
                     new Cartesian2(0, 0),
                     new Cartesian2(-1, 1),
                     new Cartesian2(-2, 2),
                     new Cartesian2(-3, 3)];

    it('create axis aligned bounding rectangle', function() {
        var rectangle = BoundingRectangle.fromPoints(positions);
        expect(rectangle.x).toEqual(-3);
        expect(rectangle.y).toEqual(-3);
        expect(rectangle.width).toEqual(6);
        expect(rectangle.height).toEqual(6);
    });

    it('fromPoints creates an empty rectangle with no positions', function() {
        var rectangle = BoundingRectangle.fromPoints();
        expect(rectangle.x).toEqual(0.0);
        expect(rectangle.y).toEqual(0.0);
        expect(rectangle.width).toEqual(0.0);
        expect(rectangle.height).toEqual(0.0);
    });

    it('fromPoints works with a result parameter', function() {
        var result = new BoundingRectangle();
        var rectangle = BoundingRectangle.fromPoints(positions, result);
        expect(rectangle).toBe(result);
        expect(rectangle.x).toEqual(-3);
        expect(rectangle.y).toEqual(-3);
        expect(rectangle.width).toEqual(6);
        expect(rectangle.height).toEqual(6);
    });

    it('fromPoints creates an empty rectangle with no positions', function() {
        var rectangle = BoundingRectangle.fromPoints();
        expect(rectangle.x).toEqual(0.0);
        expect(rectangle.y).toEqual(0.0);
        expect(rectangle.width).toEqual(0.0);
        expect(rectangle.height).toEqual(0.0);
    });

    it('fromExtent creates an empty rectangle with no extent', function() {
        var rectangle = BoundingRectangle.fromExtent();
        expect(rectangle.x).toEqual(0.0);
        expect(rectangle.y).toEqual(0.0);
        expect(rectangle.width).toEqual(0.0);
        expect(rectangle.height).toEqual(0.0);
    });

    it('create a bounding rectangle from an extent', function() {
        var extent = Extent.MAX_VALUE;
        var projection = new EquidistantCylindricalProjection(Ellipsoid.UNIT_SPHERE);
        var expected = new BoundingRectangle(extent.west, extent.south, extent.east - extent.west, extent.north - extent.south);
        expect(BoundingRectangle.fromExtent(extent, projection)).toEqual(expected);
    });

    it('fromExtent works with a result parameter', function() {
        var extent = Extent.MAX_VALUE;
        var expected = new BoundingRectangle(extent.west, extent.south, extent.east - extent.west, extent.north - extent.south);
        var projection = new EquidistantCylindricalProjection(Ellipsoid.UNIT_SPHERE);

        var result = new BoundingRectangle();
        var returnedResult = BoundingRectangle.fromExtent(extent, projection, result);
        expect(result).toBe(returnedResult);
        expect(returnedResult).toEqual(expected);
    });

    it('intersect works', function() {
        var rectangle1 = new BoundingRectangle(0, 0, 4, 4);
        var rectangle2 = new BoundingRectangle(2, 2, 4, 4);
        var rectangle3 = new BoundingRectangle(5, 5, 4, 4);
        expect(BoundingRectangle.intersect(rectangle1, rectangle2)).toEqual(Intersect.INTERSECTING);
        expect(BoundingRectangle.intersect(rectangle1, rectangle3)).toEqual(Intersect.OUTSIDE);
    });

    it('union works without a result parameter', function() {
        var rectangle1 = new BoundingRectangle(2.0, 0.0, 1.0, 1.0);
        var rectangle2 = new BoundingRectangle(-2.0, 0.0, 1.0, 2.0);
        var expected = new BoundingRectangle(-2.0, 0.0, 5.0, 2.0);
        var returnedResult = rectangle1.union(rectangle2);
        expect(returnedResult).toEqual(expected);
    });

    it('union works with a result parameter', function() {
        var rectangle1 = new BoundingRectangle(2.0, 0.0, 1.0, 1.0);
        var rectangle2 = new BoundingRectangle(-2.0, 0.0, 1.0, 2.0);
        var expected = new BoundingRectangle(-2.0, 0.0, 5.0, 2.0);
        var result = new BoundingRectangle(-1.0, -1.0, 10.0, 10.0);
        var returnedResult = rectangle1.union(rectangle2, result);
        expect(result).toBe(returnedResult);
        expect(returnedResult).toEqual(expected);
    });

    it('expand works if rectangle needs to grow right', function() {
        var rectangle = new BoundingRectangle(2.0, 0.0, 1.0, 1.0);
        var point = new Cartesian2(4.0, 0.0);
        var expected = new BoundingRectangle(2.0, 0.0, 2.0, 1.0);
        var result = rectangle.expand(point);
        expect(result).toEqual(expected);
    });

    it('expand works if rectangle needs x to grow left', function() {
        var rectangle = new BoundingRectangle(2.0, 0.0, 1.0, 1.0);
        var point = new Cartesian2(0.0, 0.0);
        var expected = new BoundingRectangle(0.0, 0.0, 3.0, 1.0);
        var result = rectangle.expand(point);
        expect(result).toEqual(expected);
    });

    it('expand works if rectangle needs to grow up', function() {
        var rectangle = new BoundingRectangle(2.0, 0.0, 1.0, 1.0);
        var point = new Cartesian2(2.0, 2.0);
        var expected = new BoundingRectangle(2.0, 0.0, 1.0, 2.0);
        var result = rectangle.expand(point);
        expect(result).toEqual(expected);
    });

    it('expand works if rectangle needs x to grow down', function() {
        var rectangle = new BoundingRectangle(2.0, 0.0, 1.0, 1.0);
        var point = new Cartesian2(2.0, -1.0);
        var expected = new BoundingRectangle(2.0, -1.0, 1.0, 2.0);
        var result = rectangle.expand(point);
        expect(result).toEqual(expected);
    });

    it('expand works if rectangle does not need to grow', function() {
        var rectangle = new BoundingRectangle(2.0, 0.0, 1.0, 1.0);
        var point = new Cartesian2(2.5, 0.6);
        var expected = new BoundingRectangle(2.0, 0.0, 1.0, 1.0);
        var result = rectangle.expand(point);
        expect(result).toEqual(expected);
    });

    it('expand works with a result parameter', function() {
        var rectangle = new BoundingRectangle(2.0, 0.0, 1.0, 1.0);
        var point = new Cartesian2(2.0, -1.0);
        var expected = new BoundingRectangle(2.0, -1.0, 1.0, 2.0);
        var result = new BoundingRectangle();
        var returnedResult = rectangle.expand(point, result);
        expect(returnedResult).toBe(returnedResult);
        expect(result).toEqual(expected);
    });

    it('static clone throws with no parameter', function() {
        expect(function() {
            BoundingRectangle.clone();
        }).toThrow();
    });

    it('static union throws with no left parameter', function() {
        var right = new BoundingRectangle(1.0, 2.0, 3.0, 4.0);
        expect(function() {
            BoundingRectangle.union(undefined, right);
        }).toThrow();
    });

    it('static union throws with no right parameter', function() {
        var left = new BoundingRectangle(1.0, 2.0, 3.0, 4.0);
        expect(function() {
            BoundingRectangle.union(left, undefined);
        }).toThrow();
    });

    it('static expand throws with no rectangle parameter', function() {
        var point = new Cartesian2();
        expect(function() {
            BoundingRectangle.expand(undefined, point);
        }).toThrow();
    });

    it('static expand throws with no point parameter', function() {
        var rectangle = new BoundingRectangle(1.0, 2.0, 3.0, 4.0);
        expect(function() {
            BoundingRectangle.expand(rectangle, undefined);
        }).toThrow();
    });

    it('static intersect throws with no left parameter', function() {
        var right = new BoundingRectangle(1.0, 2.0, 3.0, 4.0);
        expect(function() {
            BoundingRectangle.intersect(undefined, right);
        }).toThrow();
    });

    it('static intersect  throws with no right parameter', function() {
        var left = new BoundingRectangle(1.0, 2.0, 3.0, 4.0);
        expect(function() {
            BoundingRectangle.intersect(left, undefined);
        }).toThrow();
    });
});