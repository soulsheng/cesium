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

    // Bing Maps
    var bing = new Cesium.BingMapsTileProvider({
        server : 'dev.virtualearth.net',
        mapStyle : Cesium.BingMapsStyle.AERIAL,
        // Some versions of Safari support WebGL, but don't correctly implement
        // cross-origin image loading, so we need to load Bing imagery using a proxy.
        proxy : Cesium.FeatureDetection.supportsCrossOriginImagery() ? undefined : new Cesium.DefaultProxy('/proxy/')
    });

    var cb = new Cesium.CentralBody(ellipsoid);
    cb.dayTileProvider = bing;
    cb.nightImageSource = '../../Images/land_ocean_ice_lights_2048.jpg';
    cb.specularMapSource = '../../Images/earthspec1k.jpg';
    if (scene.getContext().getMaximumTextureSize() > 2048) {
        cb.cloudsMapSource = '../../Images/earthcloudmaptrans.jpg';
        cb.bumpMapSource = '../../Images/earthbump1k.jpg';
    }
    cb.showSkyAtmosphere = true;
    cb.showGroundAtmosphere = true;
    primitives.setCentralBody(cb);

    scene.getCamera().frustum.near = 1000.0;
    scene.getCamera().getControllers().addCentralBody();

    var transitioner = new Cesium.SceneTransitioner(scene, ellipsoid);

    ///////////////////////////////////////////////////////////////////////////
    // Add examples from the Sandbox here:

    /*var model = new Cesium.Model(Cesium.hellfire);
    //var model = new Cesium.Model(Cesium.facility);
    //var model = new Cesium.Model(Cesium.groundvehicle);
    model.scale = 1000.0;
    primitives.add(model);

    var camera = scene.getCamera();

    camera.position = Cesium.Cartesian3.UNIT_Z.negate();
    camera.direction = Cesium.Cartesian3.UNIT_Z;
    camera.up = Cesium.Cartesian3.UNIT_Y;
    camera.right = camera.direction.cross(camera.up);
    scene.getCamera().getControllers().get(0).spindleController.setEllipsoid(Cesium.Ellipsoid.UNIT_SPHERE);*/

    ///////////////////////////////////////////////////////////////////////////

    var time = 0;
    var amplitude = 0.00015;
    var frequency = 0.8;
    var init_amplitude = amplitude;
    var init_frequency = frequency;

    scene.setAnimation(function() {
        scene.setSunPosition(scene.getCamera().position);

        var camera = scene.getCamera();
        var up = camera.up;
        var direction = camera.direction;
        var right = camera.right;

        var yaw = amplitude * (Math.sin(2.0 * Math.PI * time * frequency) - Math.sin(6.0 * Math.PI * time * frequency) - Math.cos(2.0 * Math.PI * time * frequency)) * 1.0 / 3.0;
        var pitch = amplitude * (Math.sin(2.0 * Math.PI * time * frequency) - Math.sin(3.0 * Math.PI * time * frequency) - Math.cos(1.0 * Math.PI * time * frequency)) * 1.0 / 3.0;

        var yawrotate = Cesium.Quaternion.fromAxisAngle(up, yaw);
        var pitchrotate = Cesium.Quaternion.fromAxisAngle(right, pitch);
        var rotate = Cesium.Matrix3.fromQuaternion(yawrotate.multiply(pitchrotate));

        camera.direction = rotate.multiplyByVector(direction);
        camera.up = rotate.multiplyByVector(up);
        camera.right = camera.direction.cross(camera.up);

        time += 0.003;

        //////////////// Models
        var model1 = new Cesium.Model(Cesium.hellfire);
        var model2 = new Cesium.Model(Cesium.groundvehicle);

        if (time > 6.5 && time < 6.51) {
            var bs = new Cesium.BoundingSphere(ellipsoid.cartographicArrayToCartesianArray([Cesium.Cartographic.fromDegrees(-76.0, 35.0), Cesium.Cartographic.fromDegrees(-75.0, 35.0), Cesium.Cartographic.fromDegrees(-75.0, 36.0), Cesium.Cartographic.fromDegrees(-76.0, 36.0)]));

            model1.modelMatrix = Cesium.Matrix4.fromRotationTranslation(Cesium.Matrix3.IDENTITY, bs.center);
            model1.scale = 1000000.0;
            primitives.add(model1);

            scene.getAnimations().add({
                startValue : {
                    magnitude : 0.0
                },
                stopValue : {
                    magnitude : 12000000.0
                },
                onUpdate : function(mag) {
                    var direction = Cesium.Cartesian3.fromCartesian4(model1.modelMatrix.getColumn(3)).normalize();
                    var magnitude = mag.magnitude + ellipsoid.getMaximumRadius();
                    var center = direction.multiplyByScalar(magnitude);
                    model1.modelMatrix = Cesium.Matrix4.fromRotationTranslation(Cesium.Matrix3.IDENTITY, center);
                }
            });
        }

        if (time > 13.7 && time < 13.71) {
            var bs = new Cesium.BoundingSphere(ellipsoid.cartographicArrayToCartesianArray([Cesium.Cartographic.fromDegrees(-80.0, 30.0), Cesium.Cartographic.fromDegrees(-70.0, 30.0), Cesium.Cartographic.fromDegrees(-70.0, 40.0), Cesium.Cartographic.fromDegrees(-80.0, 40.0)]));

            model2.modelMatrix = Cesium.Matrix4.fromRotationTranslation(Cesium.Matrix3.IDENTITY, bs.center);
            model2.scale = 500000.0;
            primitives.add(model2);

            scene.getAnimations().add({
                startValue : {
                    magnitude : 0.0
                },
                stopValue : {
                    magnitude : 12000000.0
                },
                onUpdate : function(mag) {
                    var direction = Cesium.Cartesian3.fromCartesian4(model2.modelMatrix.getColumn(3)).normalize();
                    var magnitude = mag.magnitude + ellipsoid.getMaximumRadius();
                    var center = direction.multiplyByScalar(magnitude);
                    model2.modelMatrix = Cesium.Matrix4.fromRotationTranslation(Cesium.Matrix3.IDENTITY, center);
                }
            });
        }

        ///////////// Shake ramps up due to objects - all time based
        if (time > 7.0 && time < 8.0) {
            amplitude += 0.00008;
            frequency += 0.01;
        }

        if (time > 14.0 && time < 15.0) {
            amplitude += 0.00012;
            frequency += 0.012;
        }

        /////// damp frequency and amplitude
        if (amplitude > init_amplitude) {
            amplitude *= 0.98;
        }
        if (frequency > init_frequency) {
            frequency *= 0.997;
        }
    });

    (function tick() {
        scene.render();
        Cesium.requestAnimationFrame(tick);
    }());

    ///////////////////////////////////////////////////////////////////////////
    // Example mouse & keyboard handlers

    var handler = new Cesium.EventHandler(canvas);

    handler.setMouseAction(function(movement) {
        /* ... */
        // Use movement.startPosition, movement.endPosition
    }, Cesium.MouseEventType.MOVE);

    function keydownHandler(e) {
        switch (e.keyCode) {
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
        default:
            break;
        }
    }
    document.addEventListener('keydown', keydownHandler, false);

    canvas.oncontextmenu = function() {
        return false;
    };

    ///////////////////////////////////////////////////////////////////////////
    // Example resize handler

    var onResize = function() {
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