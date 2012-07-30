/*global require*/
require({
    baseUrl : '../../Source'
}, ['Cesium'], function(Cesium) {
    "use strict";
    //A real application should require only the subset of modules that
    //are actually used, instead of requiring the Cesium module, which
    //includes everything.

    var canvas = document.getElementById('glCanvas');
    var ellipsoid = Cesium.Ellipsoid.WGS84; // Used in many Sandbox examples
    var scene = new Cesium.Scene(canvas);
    var primitives = scene.getPrimitives();

    //Cesium.TerrainProvider.wireframe = true;

//    var terrainProvider = new Cesium.EllipsoidTerrainProvider(new Cesium.WebMercatorTilingScheme({
//        ellipsoid : ellipsoid,
//        numberOfLevelZeroTilesX : 2,
//        numberOfLevelZeroTilesY : 2
//    }));

    var terrainProvider = new Cesium.ArcGisImageServerTerrainProvider({
        url : 'http://elevation.arcgisonline.com/ArcGIS/rest/services/WorldElevation/DTMEllipsoidal/ImageServer',
        token : 'b_Fuz3aqYM9W6T_FcJEtgqeTsPx-hKKahTeu82QK5RWLcBZFwo-M7vAvJB4Bgkjp9by2_B1SLAXfwIdyAA-dtA..',
        proxy : new Cesium.DefaultProxy('/terrain/')
    });

    var tilingScheme = new Cesium.WebMercatorTilingScheme();
    var x = 22370;
    var y = 51108;
    var tile = new Cesium.Tile({x: x, y: y, level: 17, tilingScheme: tilingScheme});
    Cesium.when(terrainProvider.requestTileGeometry(tile), function() {
        var image = tile.geometry;
        var pixels = Cesium.getImagePixels(image);
        var vertices = new Float32Array(image.width * image.height * 3);

        var southwest = tilingScheme.cartographicToWebMercator(tile.extent.west, tile.extent.south);
        var northeast = tilingScheme.cartographicToWebMercator(tile.extent.east, tile.extent.north);
        var webMercatorExtent = {
                west : southwest.x,
                south : southwest.y,
                east : northeast.x,
                north : northeast.y
        };

        Cesium.HeightmapTessellator.computeVertices({
            heightmap: pixels,
            heightScale: 1000.0,
            heightOffset: 1000.0,
            bytesPerHeight: 3,
            strideBytes: 4,
            width: image.width,
            height: image.height,
            extent: webMercatorExtent,
            generateTextureCoordinates: false,
            relativeToCenter: new Cesium.Cartesian3(0.0, 0.0, 0.0),
            vertices: vertices,
            radiiSquared: ellipsoid.getRadiiSquared(),
            oneOverCentralBodySemimajorAxis: ellipsoid.getOneOverRadii().x
        });
        var northwestCornerX = image.width * x;
        var northwestCornerY = image.height * y;
        for (var level = 16; level >= 0; --level) {
            x /= 2;
            y /= 2;
            tile = new Cesium.Tile({x: x | 0, y: y | 0, level: level, tilingScheme: tilingScheme});
            (function(tile, x, y, level) {
                Cesium.when(terrainProvider.requestTileGeometry(tile), function() {
                    if (tile.state === Cesium.TileState.UNLOADED) {
                        return;
                    }

                    var parentImage = tile.geometry;
                    var parentPixels = Cesium.getImagePixels(parentImage);
                    var parentVertices = new Float32Array(parentImage.width * parentImage.height * 3);

                    southwest = tilingScheme.cartographicToWebMercator(tile.extent.west, tile.extent.south);
                    northeast = tilingScheme.cartographicToWebMercator(tile.extent.east, tile.extent.north);
                    webMercatorExtent = {
                            west : southwest.x,
                            south : southwest.y,
                            east : northeast.x,
                            north : northeast.y
                    };

                    Cesium.HeightmapTessellator.computeVertices({
                        heightmap: parentPixels,
                        heightScale: 1000.0,
                        heightOffset: 1000.0,
                        bytesPerHeight: 3,
                        strideBytes: 4,
                        width: parentImage.width,
                        height: parentImage.height,
                        extent: webMercatorExtent,
                        generateTextureCoordinates: false,
                        relativeToCenter: new Cesium.Cartesian3(0.0, 0.0, 0.0),
                        vertices: parentVertices,
                        radiiSquared: ellipsoid.getRadiiSquared(),
                        oneOverCentralBodySemimajorAxis: ellipsoid.getOneOverRadii().x
                    });

                    var parentNorthwestCornerX = parentImage.width * (x | 0);
                    var parentNorthwestCornerY = parentImage.height * (y | 0);

                    var maxDifference = 0.0;
                    var minDifference = 1e30;
                    var sum = 0.0;
                    var sumOfSquares = 0.0;

                    for (var v = 0; v < image.height; ++v) {
                        for (var h = 0; h < image.width; ++h) {
                            var positionIndex = (v * image.width + h) * 3;
                            var position = new Cesium.Cartesian3(vertices[positionIndex], vertices[positionIndex + 1], vertices[positionIndex + 2]);

                            var parentVFraction = (northwestCornerY + v) / (1 << (17 - level)) - parentNorthwestCornerY;
                            var parentV = parentVFraction | 0;
                            parentVFraction -= parentV;

                            if (parentV >= parentImage.height) throw new Cesium.DeveloperError();

                            var parentHFraction = (northwestCornerX + h) / (1 << (17 - level)) - parentNorthwestCornerX;
                            var parentH = parentHFraction | 0;
                            parentHFraction -= parentH;

                            if (parentH >= parentImage.width) throw new Cesium.DeveloperError();

                            var parentPositionIndex = (parentV * parentImage.width + parentH) * 3;
                            var parentPositionNW = new Cesium.Cartesian3(parentVertices[parentPositionIndex], parentVertices[parentPositionIndex + 1], parentVertices[parentPositionIndex + 2]);
                            var parentPositionNE = new Cesium.Cartesian3(parentVertices[parentPositionIndex + 3], parentVertices[parentPositionIndex + 4], parentVertices[parentPositionIndex + 5]);

                            parentPositionIndex += parentImage.width * 3;
                            var parentPositionSW = new Cesium.Cartesian3(parentVertices[parentPositionIndex], parentVertices[parentPositionIndex + 1], parentVertices[parentPositionIndex + 2]);
                            var parentPositionSE = new Cesium.Cartesian3(parentVertices[parentPositionIndex + 3], parentVertices[parentPositionIndex + 4], parentVertices[parentPositionIndex + 5]);

                            var parentPositionNorth = parentPositionNW.multiplyByScalar(1.0 - parentHFraction).add(parentPositionNE.multiplyByScalar(parentHFraction));
                            var parentPositionSouth = parentPositionSW.multiplyByScalar(1.0 - parentHFraction).add(parentPositionSE.multiplyByScalar(parentHFraction));

                            var parentPosition = parentPositionNorth.multiplyByScalar(1.0 - parentVFraction).add(parentPositionSouth.multiplyByScalar(parentVFraction));

                            var difference = parentPosition.subtract(position).magnitude();
                            sum += difference;
                            sumOfSquares += difference * difference;
                            maxDifference = Math.max(difference, maxDifference);
                            minDifference = Math.min(difference, minDifference);
                        }
                    }

                    var average = sum / (image.width * image.height);
                    var rms = Math.sqrt(sumOfSquares / (image.width * image.height));
                    console.log('Level: ' + level + ' Max error: ' + maxDifference + ' Min error ' + minDifference + " Average error: " + average + " RMS error: " + rms);
                });
            })(tile, x, y, level);
        }
    });

    var imageryLayerCollection = new Cesium.ImageryLayerCollection();

//    var esriLayer = imageryLayerCollection.addImageryProvider(new Cesium.ArcGisMapServerImageryProvider({
//        url : 'http://services.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer',
//        proxy : new Cesium.DefaultProxy('/proxy/')
//    }));

//    var streetsLayer = imageryLayerCollection.addImageryProvider(new Cesium.ArcGisMapServerImageryProvider({
//        url : 'http://server.arcgisonline.com/ArcGIS/rest/services/Reference/World_Transportation/MapServer',
//        proxy : new Cesium.DefaultProxy('/proxy/')
//    }));

    var bingAerialLayer = imageryLayerCollection.addImageryProvider(new Cesium.BingMapsImageryProvider({
        server : 'dev.virtualearth.net',
        mapStyle : Cesium.BingMapsStyle.AERIAL,
        // Some versions of Safari support WebGL, but don't correctly implement
        // cross-origin image loading, so we need to load Bing imagery using a proxy.
        proxy : Cesium.FeatureDetection.supportsCrossOriginImagery() ? undefined : new Cesium.DefaultProxy('/proxy/')
    }));

//    var bingRoadLayer = imageryLayerCollection.addImageryProvider(new Cesium.BingMapsImageryProvider({
//        server : 'dev.virtualearth.net',
//        mapStyle : Cesium.BingMapsStyle.ROAD,
//        // Some versions of Safari support WebGL, but don't correctly implement
//        // cross-origin image loading, so we need to load Bing imagery using a proxy.
//        proxy : Cesium.FeatureDetection.supportsCrossOriginImagery() ? undefined : new Cesium.DefaultProxy('/proxy/')
//    }));
//
//    var solidColorLayer = imageryLayerCollection.addImageryProvider(new Cesium.SolidColorImageryProvider());

    var testLayer = imageryLayerCollection.addImageryProvider(
            new Cesium.SingleTileImageryProvider('../../Images/TestLayer.png',
                                                 new Cesium.Extent(Cesium.Math.toRadians(-120),
                                                                   Cesium.Math.toRadians(37),
                                                                   Cesium.Math.toRadians(-119),
                                                                   Cesium.Math.toRadians(38))));

    var cb = new Cesium.CentralBody(ellipsoid, terrainProvider, imageryLayerCollection);

    cb.nightImageSource = '../../Images/land_ocean_ice_lights_2048.jpg';
    cb.specularMapSource = '../../Images/earthspec1k.jpg';
    if (scene.getContext().getMaximumTextureSize() > 2048) {
        cb.cloudsMapSource = '../../Images/earthcloudmaptrans.jpg';
        cb.bumpMapSource = '../../Images/earthbump1k.jpg';
    }
    cb.showSkyAtmosphere = true;
    cb.showGroundAtmosphere = true;
    cb.showNight = false;
    cb.affectedByLighting = false;
    primitives.setCentralBody(cb);

    scene.getCamera().frustum.near = 100.0;
    scene.getCamera().frustum.far = 20000000.0;
    scene.getCamera().getControllers().addCentralBody();

    var transitioner = new Cesium.SceneTransitioner(scene, ellipsoid);

    ///////////////////////////////////////////////////////////////////////////
    // Add examples from the Sandbox here:

    ///////////////////////////////////////////////////////////////////////////

    scene.setAnimation(function() {
        //scene.setSunPosition(scene.getCamera().position);
        scene.setSunPosition(Cesium.SunPosition.compute().position);

        // Add code here to update primitives based on changes to animation time, camera parameters, etc.
    });

    (function tick() {
        scene.render();
        Cesium.requestAnimationFrame(tick);
    }());

    ///////////////////////////////////////////////////////////////////////////
    // Example mouse & keyboard handlers

    var handler = new Cesium.EventHandler(canvas);

    var mousePosition;
    handler.setMouseAction(function(movement) {
        mousePosition = movement.endPosition;
    }, Cesium.MouseEventType.MOVE);

    (function() {
        var alphaHandler = new Cesium.EventHandler(canvas);

        var isDown = false;
        var startPosition;
        alphaHandler.setMouseAction(function(movement) {
            isDown = true;
            startPosition = movement.position;
        }, Cesium.MouseEventType.LEFT_DOWN, Cesium.EventModifier.CTRL);

        alphaHandler.setMouseAction(function(movement) {
            isDown = false;
        }, Cesium.MouseEventType.LEFT_UP, Cesium.EventModifier.CTRL);

        alphaHandler.setMouseAction(function(movement) {
            if (isDown) {
                var distance = startPosition.y - movement.endPosition.y;
                var adjustment = distance / 400;
                bingAerialLayer.alpha = Math.min(1.0, Math.max(0.0, bingAerialLayer.alpha + adjustment));
            }
        }, Cesium.MouseEventType.MOVE, Cesium.EventModifier.CTRL);
    })();

    function ellipsoidIntersections(ellipsoid, ray) {
        var position = ray.origin;
        var direction = ray.direction;
        var oneOverRadii = ellipsoid.getOneOverRadii();

        var scaledPosition = position.multiplyComponents(oneOverRadii);
        var scaledDirection = direction.multiplyComponents(oneOverRadii);

        var a = scaledDirection.magnitudeSquared();
        var b = 2.0 * scaledPosition.dot(scaledDirection);
        var c = scaledPosition.magnitudeSquared() - 1.0;

        var b2 = b * b;
        var four_ac = 4.0 * a * c;
        var radicand = b2 - four_ac;

        if (radicand < 0.0) {
            return undefined;
        }

        var q = -0.5 * (b + sign(b) * Math.sqrt(radicand));
        if (b > 0.0) {
            return [q / a, c / q];
        }
        return [c / q, q / a];
    }

    function sign(x) {
        if (x < 0.0) {
            return -1.0;
        } else if (x > 0.0) {
            return 1.0;
        }
        return 0.0;
    }

    function keydownHandler(e) {
        switch (e.keyCode) {
        case 'Q'.charCodeAt(0):
            imageryLayerCollection.raise(bingAerialLayer);
            break;
        case 'A'.charCodeAt(0):
            imageryLayerCollection.lower(bingAerialLayer);
            break;
        case 'W'.charCodeAt(0):
            imageryLayerCollection.raise(bingRoadLayer);
            break;
        case 'S'.charCodeAt(0):
            imageryLayerCollection.lower(bingRoadLayer);
            break;
        case 'E'.charCodeAt(0):
            imageryLayerCollection.raise(esriLayer);
            break;
        case 'D'.charCodeAt(0):
            imageryLayerCollection.lower(esriLayer);
            break;
        case 109: // numpad -
            imageryLayerCollection.remove(testLayer, false);
            break;
        case 107: // numpad +
            if (!imageryLayerCollection.contains(testLayer)) {
                imageryLayerCollection.add(testLayer);
            }
            break;
        case "3".charCodeAt(0): // "3" -> 3D globe
            cb.showSkyAtmosphere = true;
            cb.showGroundAtmosphere = true;
            transitioner.morphTo3D();
            break;
        case "2".charCodeAt(0): // "2" -> Columbus View
            cb.showSkyAtmosphere = false;
            cb.showGroundAtmosphere = false;
            transitioner.morphToColumbusView();
            break;
        case "1".charCodeAt(0): // "1" -> 2D map
            cb.showSkyAtmosphere = false;
            cb.showGroundAtmosphere = false;
            transitioner.morphTo2D();
            break;
        case 'F'.charCodeAt(0):
            cb._surface.toggleLodUpdate();
            break;
        case 'B'.charCodeAt(0):
            var camera = scene.getCamera();
            var pickRay = camera._getPickRayPerspective(mousePosition);
            var intersectionDistance = ellipsoidIntersections(ellipsoid, pickRay)[0];
            var intersectionPoint = pickRay.getPoint(intersectionDistance);
            var cartographicPick = ellipsoid.cartesianToCartographic(intersectionPoint);
            cb._surface.showBoundingSphereOfTileAt(cartographicPick);
            break;
        case 'V'.charCodeAt(0):
            var camera = scene.getCamera();
            console.log('position: ' + camera.getPositionWC() + ' direction: ' + camera.getDirectionWC() + ' up: ' + camera.getUpWC());
            break;
        case 'L'.charCodeAt(0):
            var position = new Cesium.Cartesian3(-2444822.7410362503, -4495391.442027786, 3800016.0770029235);
            var direction = new Cesium.Cartesian3(0.2306406390807245, 0.8212545213405562, 0.5218677100397503);
            var up = new Cesium.Cartesian3(-0.31881470184216937, -0.44294118336776583, 0.8379500545772692);
            var camera = scene.getCamera();
            camera.lookAt(position, position.add(direction), up);
            break;
        }
    }

    document.addEventListener('keydown', keydownHandler, false);

    canvas.oncontextmenu = function() {
        return false;
    };

    ///////////////////////////////////////////////////////////////////////////
    // Example resize handler

    var onResize = function () {
        var width = canvas.clientWidth;
        var height = canvas.clientHeight;

        if (canvas.width === width && canvas.height === height) {
            return;
        }

        canvas.width = width;
        canvas.height = height;

        scene.getContext().setViewport({
            x : 0,
            y : 0,
            width : width,
            height : height
        });

        scene.getCamera().frustum.aspectRatio = width / height;
    };
    window.addEventListener('resize', onResize, false);
    onResize();
});